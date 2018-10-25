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


samples = Channel.from([id: "S1221_pl",
    assembly: "${workflow.projectDir}/data/testAssembly.fasta",
    lr: "${workflow.projectDir}/data/testReads_nanopore.fastq"]
    ).view()

window = 50


// Duplicate channel
samples.into{samples_map; samples_rgi; samples_gc; samples_split}

// Take assembly and split into main chromosome and supposed plasmids
process map_longreads {
    publishDir "${params.outFolder}/${id}", mode: 'copy'

    input:
    set id, assembly, lr from samples_map


    output:
    set id, assembly, file("${id}_lr.bam"), file("${id}_lr.bai") into mos_depth

    script:
    """
    minimap2 -ax map-ont -t ${params.cpu} ${assembly} ${lr} \
    | samtools sort | samtools view -b -F 4 -o  ${id}_lr.bam 
    samtools index ${id}_lr.bam ${id}_lr.bai
    """
}


process identify_resistance_genes {
    publishDir "${params.outFolder}/${id}/rgi/", mode: 'copy'
    
    input:
    set id, assembly, lr from samples_rgi
    

    output:
    set id, file("${id}_rgi.gff3") into rgi_gff

    script:
    """
    ${RGI} -i ${assembly} -n ${params.cpu} -o ${id}_rgi

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
    Rscript 02_create_rgi_circos.R ${gff}
    """
}

process rename_annotations {
// Fixes contig names in the gff file. 
    //publishDir "${params.outFolder}/${id}/rgi/", mode: 'copy'

    input:
    set id, file(gff) from rgi_gff
    
    output:
    set id, file("${id}_rgi_fixed.gff3") into rgi_gff_fixed

    script:
    """
    sed 's/_[0-9]*//' ${gff} > ${id}_rgi_fixed.gff3
    """
}

process mos_depth{
// Calculate coverage depth
    publishDir "${params.outFolder}/${id}/depth", mode: 'copy'

    input:
    set id, assembly, aln_lr, aln_lr_idx from mos_depth

    output:
    file("${id}*")
    set id, file("${id}.regions.bed.gz") into cov_regions

    script:
    """
    source activate mosdepth
    mosdepth -t ${params.cpu} -n -b ${window} ${id} ${aln_lr} 
    """

}

process format_data_cov{
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

process calcGC{
// Calculate gc conten
    publishDir "${params.outFolder}/${id}/gc", mode: 'copy'

    input:
    set id, assembly, aln_lr, aln_lr_idx from samples_gc
    
    output:
    set id, file('gc50.txt'), file('gc1000.txt') into circos_data_gc

    script:
    """
    Rscript 01_calculate_GC.R ${assembly} 
    """
}

//Split fasta into single contigs
contigs = Channel.create()

