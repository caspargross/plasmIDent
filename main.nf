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

// Input for each sample/assembly:
// 1) ID (ESBL1991_unicycler)
// 2) Assembly (ESBL1991_final_assembly.fasta)
// 3) Long reads (ESBL1991_nanopore.fastq)


/* 

Processes:

- Map longreads with minimap2
- Find resistance genes with rgi
- Calculate coverage depth with mosdepth
- Calculate GC content with own R-Script
- Split into contigs
- Filter contigs for max/min length
- Identify gapcloser reads with own tools
- Calculate barstacking
- Write circos data
- Create circos plot

*/

ovlp = 3000
params.maxLength = 500000
params.seqPadding = 3000
window = 50

samples = Channel.from([id: "03",
    assembly: "${workflow.projectDir}/data/VRE_assembly.fasta",
    lr: "${workflow.projectDir}/data/VRE_reads.fastq"] )
    .view()

// Duplicate channel
samples.into{samples_rgi; samples_gc; samples_split}

// Split into contigs and filter for length channel
samples_split
    .map{[
        it['id'],
        file(it.get('assembly')),
        it['lr']
        ]}
    .splitFasta(record: [id: true, seqString: true])
    .map{
        def id = it[0]
        def lr = it[2]
        def contigName = it[1]['id']
        def length = it[1]['seqString'].length()
        def sequence = it[1]['seqString']
        [id, lr, contigName, length, sequence]
       }
    .filter{it[3] < params.maxLength}
    .into{contigs; contigs_2}


process pad_plasmids {
// Add prefix and suffix with sequence from oppsig end to each plasmid
    tag{id}

    input: 
    set id, lr, contigName, length, sequence from contigs_2

    output: 
    set id, file("${id}_${contigName}_padded.fasta"), lr, contigName into contigs_padded
    
    shell:
    '''
    echo '>!{contigName}' >  !{id}_!{contigName}_padded.fasta

    echo !{sequence} | awk '{print \
        substr($1, length($1)-(!{params.seqPadding} - 1), length($1))\
        $1 \
        substr($1, 1, !{params.seqPadding})\
        }' >> !{id}_!{contigName}_padded.fasta

    '''
}


// Take assembly and split into main chromosome and supposed plasmids
process map_longreads {
    publishDir "${params.outFolder}/${id}", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr, contigName from contigs_padded.groupTuple()

    output:
    set id, file("${id}_padded.fasta"), file("${id}_lr.bam"), file("${id}_lr.bai") into bam_lr

    script:
    """
    cat \$(echo ${assembly} | tr -d '[],') > ${id}_padded.fasta 
    LR=\$(echo ${lr} | tr -d '[],' | xargs -n 1 | uniq )
    minimap2 -ax map-ont -t ${params.cpu} ${id}_padded.fasta \$LR \
    | samtools sort | samtools view -b -F 4 -o  ${id}_lr.bam 
    samtools index ${id}_lr.bam ${id}_lr.bai
    """
}

bam_lr.into{bam_cov; bam_ovlp}

process find_ovlp_reads {
// Creates circos file from bam, uses R script to find overlapping reads
    tag{contig_name}

    input:
    set id, lr, contig_name, length, seq, assembly, bam, bai from contigs.combine(bam_ovlp, by : 0)
    output:
    set id, contig_name, length, file("reads.txt"), file("ovlp.txt") into circos_reads 

    script:
    """
    bedtools bamtobed -i ${bam} > reads.bed
    echo -e ${contig_name}'\\t'\$(expr ${params.seqPadding} - 10)'\\t'\$(expr ${params.seqPadding} + 10) > breaks.bed
    echo -e ${contig_name}'\\t'\$(expr ${length} - ${params.seqPadding} - 10)'\\t'\$(expr ${length} - ${params.seqPadding} + 10) >> breaks.bed
    bedtools intersect -a reads.bed -b breaks.bed -wa > ovlp.bed
    03_prepare_reads.R ovlp.bed ${contig_name} ${length} ${params.seqPadding} ovlp.txt TRUE
    03_prepare_reads.R reads.bed ${contig_name} ${length} ${params.seqPadding} reads.txt FALSE
    """
}

process identify_resistance_genes {
    publishDir "${params.outFolder}/${id}/rgi/", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr from samples_rgi
    
    output:
    set id, file("${id}_rgi.txt") into rgi_txt

    script:
    """
    source activate rgi
    rgi main -i ${assembly} -n ${params.cpu} -o ${id}_rgi
    """
}

process rename_annotations {
// Fixes contig names in the gff file. 
    input:
    set id, file(gff) from rgi_txt
    
    output:
    set id, file("${id}_rgi_fixed.txt") into rgi_gff_fixed

    script:
    """
    sed 's/_[0-9]*//' ${gff} > ${id}_rgi_fixed.txt
    """
}

process format_data_rgi {
// Converts gff file to circos readable format    
    input:
    set id, gff from rgi_gff_fixed

    output:
    set id, file("rgi.txt"), file("rgi_span.txt") into circos_data_rgi

    script:
    """
    02_create_rgi_circos.R ${gff}
    """
}

process mos_depth {
// Calculate coverage depth
    publishDir "${params.outFolder}/${id}/depth", mode: 'copy'

    input:
    set id, assembly, aln_lr, aln_lr_idx from bam_cov

    output:
    file("${id}*")
    set id, file("${id}.regions.bed.gz") into cov_regions

    script:
    """
    source activate mosdepth
    mosdepth -t ${params.cpu} -n -b ${window} ${id} ${aln_lr} 
    """

}

process format_data_cov {
// Formats coverage data for use in circos
    
    input:
    set id, bed from cov_regions

    output:
    set id, file("cov.bed") into circos_data_cov

    script:
    """
    gunzip -c ${bed} > cov.bed
    """
}

process calcGC {
// Calculate gc conten
   publishDir "${params.outFolder}/${id}/gc", mode: 'copy'

    input:
    set id, assembly, lr from samples_gc
    
    output:
    set id, file('gc50.txt'), file('gc1000.txt') into circos_data_gc

    script:
    """
    01_calculate_GC.R ${assembly} 
    """
}

// Combine all finished circos data based on the id
circos_data_gc
   .join(circos_data_cov)
       .join(circos_data_rgi)
       .set{circos_data}
// id | gc50 | gc1000 | cov | rgi | rgi_span


// Combine contig data with sample wide circos data
combined_data = circos_reads.combine(circos_data, by: 0)

process circos{
// Use the combined data to create nice circos plots
publishDir "${params.outFolder}/${id}/plasmidPlots", mode: 'copy'

input:
set id, contigID, length, file(reads), file(ovlp),  file(gc50), file(gc1000), file(cov), file(rgi), file(rgi_span) from combined_data

output:
file("${id}_${contigID}_plasmid.*")

script:
"""
echo "chr	-	${contigID}	1	0	${length}	chr1	color=lblue" > contig.txt
ln -s ${workflow.projectDir}/circos_confs//* .
circos
mv circos.png ${id}_${contigID}_plasmid.png
mv circos.svg ${id}_${contigID}_plasmid.svg
"""
}

