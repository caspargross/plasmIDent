/**
*
*   RESITANCE PLASMID IDENTIFICATION PIPELINE
*
*   Nextflow pipeline to analyse assembled bacterial genomes 
*   using long reads to identify circular plasmids with and
*   antibiotical resistance genes
*   Caspar Gross 2018
* 
**/

// Run parameter:
params.assembler = "unicycler"
params.folder_suffix = "_all"

// Arguments:
files = Channel.fromPath(params.pathfile)
    .ifEmpty {error "Cannot find file with path locations in ${params.pathfile}"}
    .splitCsv(header: true)
    .view()


// Filter contigs by length:
process contig_length_table {
    tag{id}
    publishDir "${params.outFolder}/${id}/"

    input:
    set id, sr1, sr2, lr from files
    
    output:
    set id, sr1, sr2, lr, file("plasmid*.fasta") into plasmids

    script:
    """
    #!/usr/bin/Rscript

    library(seqinr)
    library(data.table)

    genome <- read.fasta("${params.inFolder}/${id}${params.folder_suffix}/${id}_${params.assembler}_final.fasta")
    
    # Calculate contig lengths
    dt <- data.table(file='${id}_${params.assembler}',
        length = lapply(genome, length),
        contig = lapply(genome, function(x) attributes(x)\$name)
    )

    # Write contig lengths to file
    fwrite(dt, file = 'contig_lengths.csv')
     
    # Filter contigs for length and extract into new file
    selection <- dt[length>${params.minLength} & length<${params.maxLength},,]

    for (i in 1:nrow(selection)) {
        contig = genome[selection\$contig[[i]]]
        write.fasta(sequences=contig, names=names(contig), file.out=paste0('plasmid',i,'.fasta'))
    }

    """

}

// Duplicate channel
plasmids.into{plasmids_lr; plasmids_sr}

// Take assembly and split into main chromosome and supposed plasmids
process map_longreads {
    publishDir "${params.outFolder}/${id}"

    input:
    set id, sr1, sr2, lr, contig from plasmids_lr.transpose()

    file_name = "aln_" + contig.baseName + "_" + id + "_lr"

    output:
    file(file_name + '.bam') into aln_lr
    file(file_name + '.bai')

    script:
    """
    ${MINIMAP2} -ax map-ont -t ${params.cpu} ${contig} ${lr} \
    | ${SAMTOOLS} sort  -o ${file_name}.bam 
    ${SAMTOOLS} index ${file_name}.bam ${file_name}.bai
    """
}

process map_shortreads {
    publishDir "${params.outFolder}/sr_alignment/"

    input:
    set id, sr1, sr2, lr, plasmid from contigs_sr.transpose()

    output:
    file('alignment_lr.bam') into mappedShortReads
    file('alignment_lr.bai')

    script:
    """
    ${BWA} index ${plasmid}
    ${BWA} aln ${plasmid} ${sr1} > R1.sai
    ${BWA} aln ${plasmid} ${sr2} > R2.sai
    ${BWA} sampe ${plasmid} R1.sai R2.sai ${sr1} ${sr2} \
    | ${SAMTOOLS} sort -o alignment_sr.bam
    ${SAMTOOLS} index alignment_sr.bam alignment_sr.bai
    """

}

process identify_resistance_genes {
    publishDir "${params.outFolder}/rgi/"
    
    input:
    

    output:
    file("resistances.gff3") into rgi_gff

    script:
    """
    ${RGI} -i ${params.genome} -n ${params.cpu} -o resistances

    """
}

// The gff file produced by rgi has wrong contig names. 
// After this script they can be load directly into a genome viewer.
process rename_annotations {
    publishDir "${params.outFolder}/rgi/"

    input:
    file(gff) from rgi_gff
    
    output:
    file("res_annot.gff3")

    script:
    """
    sed 's/_[0-9]*//' ${gff} > res_annot.gff3
    """
}
