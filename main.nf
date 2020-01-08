#!/usr/bin/env nextflow

/*
Nextflow Bismark Alignment
----------------------------------------
*/

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Check Sample sheet validity
checkCSV(params.sampleSheet, 6)

// Check and import fastq files
fq = getFastq(params.sampleSheet) // Copy fastq

// Copy files to working directory - this is in $HOME
// This is not a desireable way of doing this but it works...
process copyFiles {
    tag { filename + ' - Copy files' }

    publishDir "${copy_path}", mode: 'copy'

    input:
    set group,
        sample,
        filename,
        copy_path,
        file(reads) from fq

    output:
    set group,
        sample,
        filename,
        copy_path,
        file("${filename}_R1_001.fastq.gz"),
        file("${filename}_R2_001.fastq.gz") optional true into fq_channel

    script:
    """
    echo "copying files"
    """
}

// Split fastq files into separate input channels
( ch_fastqc, ch_trimgalore ) = fq_channel.into(2)

process runFastqc {
    tag { filename + ' - Fastqc' }

    publishDir "${params.outDir}/${filename}/fastqc", mode: 'copy'

    input:
    set group,
        sample,
        filename,
        copy_path,
        file(r1),
        file(r2) from ch_fastqc

    output:
    file "*.{html,zip}" into results_fastqc

    script:
    """
    fastqc \
        -t ${task.cpus} \
        -o . \
        ${r1} ${r2}
    """
}

process runTrimgalore {

    tag { filename + ' - Trimgalore' }

    publishDir "${params.outDir}/${filename}/trimgalore", mode: 'copy'

    input:
    set group,
        sample,
        filename,
        copy_path,
        file(r1),
        file(r2) from ch_trimgalore

    output:
    file "${filename}*.{html,zip,gz,txt}" into results_trimgalore
    set filename,
        file("${filename}_R1_001_val_1.fq.gz"),
        file("${filename}_R2_001_val_2.fq.gz") optional true into trimmed_channel

    script:
    """
    /homes/nader.aryamanesh/apps_nader/TrimGalore-master/trim_galore \
        --cores ${task.cpus} \
        --path_to_cutadapt cutadapt \
        -o . \
        -paired --clip_R1 6 --clip_R2 6 -a CACTCACCAATCT -a2 CACTCACCAATCT --fastqc \
        ${r1} ${r2}
    """
}

( ch_BismarkPE ) = trimmed_channel.into(1)


process runBismarkPE {
    tag { filename + ' - BismarkPE' }

    publishDir "${params.outDir}/${filename}/bismark", mode: 'copy'

    input:
    set filename,
        file(val_1),
        file(val_2) from ch_BismarkPE

    output:
    file "${filename}*.{bam,gz,txt,log}" into results_BismarkPE
    set filename,
        file("${filename}_R1_001_val_1.fq.gz_unmapped_reads_1.fq.gz") optional true into bismarkR1_channel
    set filename,
        file("${filename}_R2_001_val_2.fq.gz_unmapped_reads_2.fq.gz") optional true into bismarkR2_channel
    set filename,
        file("${filename}_R1_001_val_1_bismark_bt2_pe.bam") optional true into bamPE_channel

    script:
    """
    echo -e "Part 1 - Running bismark_alignment for PE directional on both R1 and R2:\n"

    /homes/nader.aryamanesh/apps_nader/Bismark-master/bismark \
        -N 1 -L 20 --un \
        --parallel 2 \
        -p ${task.cpus} \
        /homes/nader.aryamanesh/projects/191111_pig_embryo_wgbs/reference \
        -1 ${val_1} -2 ${val_2}
    """
}

process runBismarkR1 {
    tag { filename + ' - BismarkR1' }

    publishDir "${params.outDir}/${filename}/bismark", mode: 'copy'

    input:
    set filename,
        file(un_fq_1) from bismarkR1_channel

    output:
    file "${filename}*.{bam,gz,txt,log}" into results_BismarkR1
    set filename,
        file("${filename}_R1_001_val_1.fq.gz_unmapped_reads_1_bismark_bt2.bam") optional true into bamR1_channel

    script:
    """
    echo -e "Part 2 - Running bismark_alignment for SE directional on unmapped R1:\n"

    /homes/nader.aryamanesh/apps_nader/Bismark-master/bismark \
    -N 1 -L 20 \
    --parallel 2 \
    -p ${task.cpus} \
    /homes/nader.aryamanesh/projects/191111_pig_embryo_wgbs/reference \
    ${un_fq_1}
    """
}

process runBismarkR2 {
    tag { filename + ' - BismarkR2' }

    publishDir "${params.outDir}/${filename}/bismark", mode: 'copy'

    input:
    set filename,
        file(un_fq_2) from bismarkR2_channel

    output:
    file "${filename}*.{bam,gz,txt,log}" into results_BismarkR2
    set filename,
        file("${filename}_R2_001_val_2.fq.gz_unmapped_reads_2_bismark_bt2.bam") optional true into bamR2_channel

    script:
    """
    echo -e "Part 3 - Running bismark_alignment for SE directional pbat on unmapped R2:\n"

    /homes/nader.aryamanesh/apps_nader/Bismark-master/bismark \
        -N 1 -L 20 --pbat \
        --parallel 2 \
        -p ${task.cpus} \
        /homes/nader.aryamanesh/projects/191111_pig_embryo_wgbs/reference \
        ${un_fq_2}
    """
}

( ch_Bismark_extractor_bamPE, ch_Samtools_bamPE ) = bamPE_channel.into(2)
( ch_Bismark_extractor_bamR1, ch_Samtools_bamR1 ) = bamR1_channel.into(2)
( ch_Bismark_extractor_bamR2, ch_Samtools_bamR2 ) = bamR2_channel.into(2)


process runBismark_extractor {
    tag { filename + ' - BismarkExtractor' }

    publishDir "${params.outDir}/${filename}/bismarkExtractor", mode: 'copy'

    input:
    set filename,
        file(bamPE) from ch_Bismark_extractor_bamPE
    set filename,
        file(bamR1) from ch_Bismark_extractor_bamR1
    set filename,
        file(bamR2) from ch_Bismark_extractor_bamR2

    output:
    file "${filename}*.{bam,gz,txt,log}" optional true into results_BismarkExtractor

    script:
    """
		echo -e "Methylation extraction PE\n"
    /homes/nader.aryamanesh/apps_nader/Bismark-master/bismark_methylation_extractor \
        --comprehensive --report --buffer_size 8G --paired-end --no_overlap --gzip --bedGraph \
        --parallel ${task.cpus} \
        ${bamPE}

		echo -e "Methylation extraction SE1\n"
    /homes/nader.aryamanesh/apps_nader/Bismark-master/bismark_methylation_extractor \
        --comprehensive --report --buffer_size 8G -s --gzip --bedGraph \
        --parallel ${task.cpus} \
        ${bamR1}

    echo -e "Methylation extraction SE1\n"
    /homes/nader.aryamanesh/apps_nader/Bismark-master/bismark_methylation_extractor \
        --comprehensive --report --buffer_size 8G -s --gzip --bedGraph \
        --parallel ${task.cpus} \
        ${bamR2}
    zcat CpG*.txt.gz > ${filename}_CpG_context_merged.txt
    zcat CHG*.txt.gz > ${filename}_CHG_context_merged.txt
    zcat CHH*.txt.gz > ${filename}_CHH_context_merged.txt
    gzip ${filename}_CpG_context_merged.txt
    gzip ${filename}_CHG_context_merged.txt
    gzip ${filename}_CHH_context_merged.txt
    """
}


process runSamtools_view {
    tag { filename + ' - samtools_view' }

    publishDir "${params.outDir}/${filename}/samtools", mode: 'copy'

    input:
    set filename,
        file(bamPE) from ch_Samtools_bamPE
    set filename,
        file(bamR1) from ch_Samtools_bamR1
    set filename,
        file(bamR2) from ch_Samtools_bamR2

    output:
    set filename,
        file("${filename}_bismark_PE.bam"),
        file("${filename}_unmapped_bismark_SE1.bam"),
        file("${filename}_unmapped_bismark_SE2.bam") optional true into merge_channel

    script:
    """
    echo -e "Running samtools view"
		samtools view -@ ${task.cpus} -b -h ${bamPE} > ${filename}_bismark_PE.bam
		samtools view -@ ${task.cpus} -b -h ${bamR1} > ${filename}_unmapped_bismark_SE1.bam
		samtools view -@ ${task.cpus} -b -h ${bamR2} > ${filename}_unmapped_bismark_SE2.bam
    """
}


process runSamtools_merge {
    tag { filename + ' - samtools_merge' }

    publishDir "${params.outDir}/${filename}/samtools", mode: 'copy'

    input:
    set filename,
        file(bamPE),
        file(bamR1),
        file(bamR2) from merge_channel

    output:
    set filename,
        file("${filename}_bismark_combined.bam") optional true into sam_merge_channel

    script:
    """
    echo -e "combining bam files"
		samtools merge -@ ${task.cpus} -h ${bamPE} ${filename}_bismark_combined.bam ${bamPE} ${bamR1} ${bamR2}
    """
}

( ch_Samtools_merge_sort, ch_Bismark_extractor_merged ) = sam_merge_channel.into(2)

process runSamtools_merge_sort {
    tag { filename + ' - samtools_sort' }

    publishDir "${params.outDir}/${filename}/samtools", mode: 'copy'

    input:
    set filename,
        file(bam_merge) from ch_Samtools_merge_sort

    output:
    set filename,
        file("${filename}_bismark_combined_sort.bam") optional true into sam_sort_channel

    script:
    """
    samtools sort -@ ${task.cpus} -o ${filename}_bismark_combined_sort.bam ${bam_merge}
    """
}



process runSamtools_merge_index {
    tag { filename + ' - samtools_index' }

    publishDir "${params.outDir}/${filename}/samtools", mode: 'copy'

    input:
    set filename,
        file(bam_sort) from ch_Samtools_merge_index

    output:
    file "${filename}*.{bai}" into results_Samtools_index

    script:
    """
    samtools index ${bam_sort}
    """
}

process runBismark_extractor_merged {
    tag { filename + ' - BismarkExtractorMerged' }

    publishDir "${params.outDir}/${filename}/bismarkExtractorMerge", mode: 'copy'

    input:
    set filename,
        file(bamPE) from ch_Bismark_extractor_merged

    output:
    file "${filename}*.{bam,gz,txt,log}" into results_BismarkExtractor_merged

    script:
    """
		echo -e "Methylation extraction PE\n"
    /homes/nader.aryamanesh/apps_nader/Bismark-master/bismark_methylation_extractor \
        --comprehensive --report --buffer_size 8G --paired-end --no_overlap --gzip --bedGraph \
        --parallel ${task.cpus} \
        ${bamPE}

    mv CpG*.txt.gz ${filename}_CpG_context_merged.txt.gz
    mv CHG*.txt.gz ${filename}_CHG_context_merged.txt.gz
    mv CHH*.txt.gz ${filename}_CHH_context_merged.txt.gz
    """
}

/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
    println """Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """
}

/*
Pipeline functions
*/

// Check sample sheet validity - this expects a sample sheet
// following the specifications defined in the example sample
// sheet in this directory.
def checkCSV(sampleSheet, nCol) {
    // Check sample sheet exists and isn't empty
    ss = file(sampleSheet)
    if (!ss.exists()) exit 1, "Sample sheet file doesn't exist: ${sampleSheet}"
    if (ss.isEmpty()) exit 1, "Sample sheet file is empty: ${sampleSheet}"

    // Check number of columns in CSV
    Channel
        .fromPath(sampleSheet)
        .splitCsv(header: true)
        .subscribe {val ->

            // Check that the sample sheet has correct no. fields
            if (val.size() != nCol) exit 1, "CSV row is incorrect: ${val}"

            // Check CSV fields are complete (no empty fields)
            val.each { entry ->
                if(!entry.value) exit 1, "CSV value missing - Key: ${entry.key} Value: ${entry.value}"
            }
        }

    return true // Return statement

}

// Import fastq files using samplesheet
// This will build the path to the files from the
// columns in the sample sheet.
def getFastq(sampleSheet) {

    // Read lines of sampleSheet
    ssFile = file(sampleSheet)
    fqChannel = Channel
                    .from(ssFile)
                    .splitCsv(header: true)
                    .map { row ->

                        // Building file path
                        r1 = file( [ row.path, row.group, row.sample, row.R1 ].join("/") )
                        r2 = file( [ row.path, row.group, row.sample, row.R2 ].join("/") )

                        // Path
                        // pth = file([ row.path, row.group, row.sample ].join("/"))

                        // Check if reads files exist
                        if (!r1.exists()) exit 1, "Reads file doesn't exist: ${r1}"
                        if (!r2.exists()) exit 1, "Reads file doesn't exist: ${r2}"

                        // Check if file is empty
                        if (r1.isEmpty()) exit 1, "Reads file is empty: ${r1}"
                        if (r1.isEmpty()) exit 1, "Reads file is empty: ${r2}"

                        // Make copy path
                        cpy_path = file("${workflow.workDir}/fastq_copy/${row.sample}")
                        cpy_path.mkdirs()

                        // Returning list object
                        lst = [ row.group,
                                row.sample,
                                row.filename,
                                cpy_path,
                                tuple(r1, r2)
                                ]

                        return lst
                    }
}

def helpMessage() {
    message = """
    This nextflow pipeline runs the following packages:
        1- fastqc
        2- trimgalore
        3- bismark alignment
        4- bismark extractor
        5- samtools to merge the bam files

    Usage:

    The typical command for running the pipeline is as follows:

    conda activate nextflow-na
    nextflow run main_v01.nf -profile slurm -w PATH_TO_/workDir --outDir PATH_TO_/outDir

    Mandatory arguments:
      -profile                      slurm [only optimised for slurm]
      -w                            Working directory where the initial results will be saved

    Optional arguments:
      --outDir                      The output directory where the results will be saved [nf_results]
      --sampleSheet                 The PATH/NAME of the sample sheet (e.g. sampleSheet.csv)
      --userEmail                   The email address of the user[Nader.Aryamanesh@sahmri.com]
      --help                        Display this help message [False]:

      To display this help message, please type:
      nextflow run main_v03.nf --help

      For more information, contact Nader.Aryamanesh@sahmri.com
    """
    println(message)
}
