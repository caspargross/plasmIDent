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



params.maxLength = 500000
params.seqPadding = 3000
params.covWindow = 50
params.input = "plasmid_files_test.tsv"

samples = getFiles(params.input)

// Duplicate channel
samples.into{samples_rgi; samples_gc; samples_split; samples_map; samples_table}

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
  //.view()
    .into{contigs; contigs_2}


process pad_plasmids {
// Add prefix and suffix with sequence from oppsig end to each plasmid
    tag{contigName}

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

process combine_padded_contigs {
// Recombines padded contigs into a single fasta
    tag{contigName}

    input:
    set id, assembly, lr, contigName from contigs_padded.groupTuple()

    output:
    set id, file("${id}_padded.fasta"), lr, val("padded") into map_padded

    script:
    """
    cat \$(echo ${assembly} | tr -d '[],') > ${id}_padded.fasta 
    """
}

// Mix channel with padded and normal contigs
samples_map
  //.view()
    .map{[it['id'], 
        it['assembly'], 
        it['lr'], 
        'normal']}
    .mix(map_padded
        .map{[it[0], 
            it[1], 
            it[2][1],
            it[3]]})
  //.view()
    .set{to_mapping}

process map_longreads {
// Use minimap2 to align longreads to padded contigs
    publishDir "${params.outFolder}/${id}/alignment/", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr, type from to_mapping

    output:
    set id, assembly, type, file("${id}_${type}_lr.bam"), file("${id}_${type}_lr.bai") into bam_lr

    script:
    """
    minimap2 -ax map-ont -t ${params.cpu} ${assembly} ${lr} \
    | samtools sort | samtools view -b -F 4 -o  ${id}_${type}_lr.bam 
    samtools index ${id}_${type}_lr.bam ${id}_${type}_lr.bai
    """
}

// Distribute bamfiles for coverage and read overlap identification
bam_cov = Channel.create()
bam_ovlp = Channel.create()
bam_lr.into{bam_cov; bam_ovlp}

process find_ovlp_reads {
// Creates circos file from bam, uses R script to find overlapping reads
    tag{contig_name}

    input:
    set id, lr, contig_name, length, seq, assembly, type, bam, bai from contigs.combine(bam_ovlp.filter{it[2] == 'padded'}, by : 0)

    output:
    set id, contig_name, length, file("reads.txt"), file("ovlp.txt") into circos_reads 

    script:
    """
    bedtools bamtobed -i ${bam} > reads.bed
    echo -e ${contig_name}'\\t'\$(expr ${params.seqPadding} - 10)'\\t'\$(expr ${params.seqPadding} + 10) > breaks.bed
    echo -e ${contig_name}'\\t'\$(expr ${length} - ${params.seqPadding} - 10)'\\t'\$(expr ${length} - ${params.seqPadding} + 10) >> breaks.bed
    bedtools intersect -a reads.bed -b breaks.bed -wa > ovlp.bed
    03_prepare_bed.R ovlp.bed ${params.seqPadding} ovlp.txt ${contig_name} ${length} TRUE
    03_prepare_bed.R reads.bed ${params.seqPadding} reads.txt  ${contig_name} ${length} FALSE
    """
}

process identify_resistance_genes {
// Find antibiotic resistance genes in the CARD database
    publishDir "${params.outFolder}/${id}/resistances", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr from samples_rgi
    
    output:
    set id, file("${id}_rgi.txt") into from_rgi

    script:
    """
    source activate rgi
    rgi main -i ${assembly} -n ${params.cpu} -o ${id}_rgi
    """
}

from_rgi.into{rgi_txt; table_data_rgi}

process format_data_rgi {
// Converts gff file to circos readable format    
    tag{id}

    input:
    set id, rgi from rgi_txt

    output:
    set id, file("rgi.txt"), file("rgi_span.txt") into circos_data_rgi

    script:
    """
    02_create_rgi_circos.R ${rgi}
    """
}

process mos_depth {
// Calculate coverage depth
    publishDir "${params.outFolder}/${id}/coverage", mode: 'copy'
    tag{id}

    input:
    set id, assembly, type, aln_lr, aln_lr_idx from bam_cov

    output:
    file("${id}*")
    set id, file("${id}.regions.bed.gz"), type into cov_bed

    script:
    """
    source activate mosdepth
    mosdepth -t ${params.cpu} -n -b ${params.covWindow} ${id} ${aln_lr} 
    """
}

process format_data_cov {
// Formats coverage data for use in circos
    tag{id}
    
    input:
    set id, bed, type from cov_bed

    output:
    set id, file("cov.txt"), type into cov_formated

    script:
    if (type == "padded")
        """
        gunzip -c ${bed} > cov.bed
        03_prepare_bed.R cov.bed  ${params.seqPadding} cov.txt
        """
    else
        """
        gunzip -c ${bed} > cov.txt
        """
}

// Distribute coverage file for circos (padded)  and summary table (normal)
circos_data_cov = Channel.create()
table_data_cov = Channel.create()
cov_formated.choice(circos_data_cov, table_data_cov) { it[2] == 'padded' ? 0 : 1 }

process calcGC {
// Calculate gc conten
    publishDir "${params.outFolder}/${id}/gc", mode: 'copy'
    tag{id}

    input:
    set id, assembly, lr from samples_gc
    
    output:
    set id, file('gc1000.txt'), assembly into table_data_gc
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

// Combine contig data with sample wide circos data
combined_data = circos_reads.combine(circos_data, by: 0)

// Combine all table data based on id
table_data_gc
    .join(table_data_cov)
        .join(table_data_rgi)
        .set{table_data}

process circos{
// Use the combined data to create circular plots
    publishDir "${params.outFolder}/${id}/plots", mode: 'copy'
    tag{contigID}

    input:
    set id, contigID, length, file(reads), file(ovlp), file(gc50), file(gc1000), file(cov), type, file(rgi), file(rgi_span) from combined_data

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

process table{
// Create table with contig informations
    publishDir "${params.outFolder}/${id}/", mode: 'copy'

    input:
    set id, gc, assembly, cov, type, rgi from table_data

    output:
    file("${id}_summary.csv")

    script:
    """
    04_summary_table.R ${assembly} ${rgi} ${cov} ${gc}
    mv contig_summary.txt ${id}_summary.csv
    """
}


/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def getFiles(tsvFile) {
  // Extracts Read Files from TSV
  log.info "Read input file: " + tsvFile
  log.info "------------------------------"
  Channel.fromPath(tsvFile)
      .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
      .splitCsv(sep:'\t', skip: 1)
      .map { row ->
            //def id = row[0]
            //def assembly = returnFile(row[1])
            //def lr = returnFile(row[2])
            [id:row[0], assembly:returnFile(row[1]), lr:returnFile(row[2])]
            }
       .view()
}

def returnFile(it) {
// Return file if it exists
    if (workflow.profile in ['test', 'localtest'] ) {
        inputFile = file("$workflow.projectDir/" + it)
    } else {
        inputFile = file(it)
    }
    if (!file(inputFile).exists()) exit 1, "Missing file in TSV file: ${inputFile}, see --help for more information"
    return inputFile
}

