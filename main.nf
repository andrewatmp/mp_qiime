params.reads = "$projectDir/seqs"
params.fastq = "$projectDir/seqs/*.fastq.gz"
params.trunclen = 415
params.minreads = 100
params.refseqs = "$projectDir/ncbi-refseqs.qza"
params.reftax =  "$projectDir/ncbi-refseqs-taxonomy.qza"
params.maxaccepts = 1
params.artifacts = "$projectDir/artifacts"
params.outputdir = "$projectDir/multiqc"

log.info """\
    MP - Q I I M E   P I P E L I N E
    ===================================
    Reads        : ${params.reads}
    Trunc Len    : ${params.trunclen}
    Min Reads    : ${params.minreads}
    """
    .stripIndent(true)

println "reads: $params.reads"

process FASTQC {

    tag "FastQC"
    container "andrewatmp/testf"
    label "big_memory"

    input:
    path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc $reads
    """
}

process IMPORT {
    tag "Importing sequences"
    // publishDir "${params.artifacts}", mode: 'copy'
    container "andrewatmp/testf"

    input:
    path(reads)

    output:
    path("demux.qza"), emit: demux
    path("demux.qzv"), emit: demuxvis

    script:

    """
    qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path ${reads} \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path demux.qza

    qiime demux summarize \
    --i-data demux.qza \
    --o-visualization demux.qzv
    """
}

process DEMUXVIS {

    tag "Quality Visualization"
    container "andrewatmp/testf"


    input:
    path(demux)

    output:
    path("demux.qzv"), emit: demuxvis

    script:

    """
    qiime demux summarize \
    --i-data ${demux} \
    --o-visualization demuxvis.qzv
    """
}

process DADA {

    tag "Dada2 Error Correction"
    container "andrewatmp/testf"


    input:
    path(qza)
    
    output:
    path("rep-seqs.qza"), emit: repseqs
    path("table.qza"), emit: table
    path("stats.qza"), emit: stats

    script:

    """
    qiime dada2 denoise-pyro \
    --i-demultiplexed-seqs $qza \
    --p-trunc-len ${params.trunclen} \
    --p-trim-left 15 \
    --p-trunc-q 1 \
    --p-max-ee 4 \
    --o-representative-sequences rep-seqs.qza \
    --o-table table.qza \
    --o-denoising-stats stats.qza \
    --verbose
    """

}

process MINREADS {

    tag "Filtering for min reads"
    container "andrewatmp/testf"


    input:
    path(table)

    output:
    path("filtered-table.qza"), emit: filtered

    script:

    """
    qiime feature-table filter-features \
    --i-table ${table} \
    --p-min-frequency ${params.minreads} \
    --o-filtered-table filtered-table.qza
    """
}

process DADARESULTS {

    tag "Generate dada visualizations"
    container "andrewatmp/testf"


    input:
    path(repseqs)
    path(table)
    path(stats)
    path(filtered)

    output:
    path("rep-seqs.qzv"), emit: repseqsvis
    path("table.qzv"), emit: tablevis
    path("stats.qzv"), emit: statsvis
    path("filtered-table.qzv"), emit: filteredtablevis


    script:

    """
    qiime feature-table tabulate-seqs \
    --i-data $repseqs \
    --o-visualization rep-seqs.qzv

    qiime feature-table summarize \
    --i-table $table \
    --o-visualization table.qzv

    qiime feature-table summarize \
    --i-table $filtered \
    --o-visualization filtered-table.qzv

    qiime metadata tabulate \
    --m-input-file $stats \
    --o-visualization stats.qzv
    """
}

process CLASSIFY {

    tag "Classify using BLAST"
    container "andrewatmp/testf"


    input:
    path(refseqs)
    path(reftax)
    path(repseqs)

    output:
    path("classification.qza"), emit: classification
    path("blastresults.qza"), emit: blastresults

    script:

    """
    qiime feature-classifier classify-consensus-blast \
    --i-query $repseqs \
    --i-reference-reads $refseqs \
    --i-reference-taxonomy $reftax \
    --p-maxaccepts ${params.maxaccepts} \
    --p-perc-identity 0.99 \
    --o-classification classification.qza \
    --o-search-results blastresults.qza 
    """
}

process TABULATE {

    tag "Tabulate Classify Results"
    container "andrewatmp/testf"


    input:
    path(classification)
    path(blastresults)

    output:
    path("classification.qzv"), emit: classificationvis
    path("blastresults.qzv"), emit: blastresultsvis
    
    script:
    """
    qiime metadata tabulate \
    --m-input-file $blastresults \
    --o-visualization blastresults.qzv

    qiime metadata tabulate \
    --m-input-file $classification \
    --o-visualization classification.qzv
  """
}

process BARPLOT {

    tag "Generate barplot"
    container "andrewatmp/qiime_unzip"

    input:
    path(filtered)
    path(classification)

    output:
    path("taxa-bar-plots.qzv"), emit: barplot
    path("data/*"), emit: data

    script:

    """
    qiime taxa barplot \
    --i-table $filtered \
    --i-taxonomy $classification \
    --o-visualization "taxa-bar-plots.qzv"

    mkdir extracted
    mkdir data
    unzip taxa-bar-plots.qzv '*/data/*' -d extracted
    mv extracted/*/data/* data
    mv data/index.html data/Taxonomy_mqc.html
    rm -rf extracted
    """

}

process MULTIQC {

    tag "MultiQC"
    container "andrewatmp/multiqc"
    stageInMode 'copy'


    input:
    path(fastqc)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

process MULTIQC2 {

    tag "MultiQC2"
    container "andrewatmp/multiqc"

    input:
    path(data)

    script:
    """
    multiqc .
    """
}




workflow {

    fastqc_ch=Channel.fromPath(params.fastq)
    FASTQC(fastqc_ch)
    IMPORT(params.reads)

    dada_ch = DADA(IMPORT.out.demux)
    filtered_ch = MINREADS(DADA.out.table)
    DADARESULTS(dada_ch, filtered_ch)

    classification_ch = CLASSIFY(params.refseqs, params.reftax, DADA.out.repseqs)
    TABULATE(classification_ch)

    BARPLOT(filtered_ch, CLASSIFY.out.classification)

    multiqc_files = Channel.empty()
    multiqc_files = multiqc_files.mix(FASTQC.out.fastqc_results)
    multiqc_files = multiqc_files.mix(BARPLOT.out.data)
    MULTIQC(multiqc_files.collect())


    // MULTIQC2(BARPLOT.out.data)
}
