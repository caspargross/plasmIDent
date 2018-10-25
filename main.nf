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


samples = Channel.from([id: "ESBl1991",
    assembly: "${workflow.projectDir}/data/testAssembly.fasta",
    lr: "${workflow.projectDir}/data/testReads_nanopore.fastq"]
    ).view()

window = 50


// Duplicate channel
samples.into{samples_map; samples_rgi}

// Take assembly and split into main chromosome and supposed plasmids
process map_longreads {
    publishDir "${params.outFolder}/${id}"

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
    publishDir "${params.outFolder}/rgi/"
    
    input:
    set id, assembly, lr from samples_rgi
    

    output:
    file("resistances.gff3") into rgi_gff

    script:
    """
    ${RGI} -i ${assembly} -n ${params.cpu} -o resistances

    """
}

// The gff file produced by rgi has wrong contig names. 
// After this script they can be load directly into a genome viewer.
process rename_annotations {
    publishDir "${params.outFolder}/rgix/"

    input:
    file(gff) from rgi_gff
    
    output:
    file("res_annot.gff3")

    script:
    """
    sed 's/_[0-9]*//' ${gff} > res_annot.gff3
    """
}

process mos_depth{
    publishDir "${params.outFolder}/depth"

    input:
    set id, assembly, aln_lr, aln_lr_idx from  mos_depth

    output:
    file("${id}*")

    script:
    """
    source activate mosdepth
    mosdepth -t ${params.cpu} -n -b ${window} ${id} ${aln_lr} 
    """

}
