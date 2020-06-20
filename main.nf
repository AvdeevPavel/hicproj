#!/usr/bin/env nextflow

def helpMessage() {
    // Add to this help message with new command line parameters
    // loger info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --libs library.yml --graph graph.gfa --output_dir sb100

    Mandatory arguments:
        --hic_lib                          Path to library with Hi-C data (must be surrounded with quotes)
        --graph                            Path to input assembly graph in GFA format (must be surrounded with quotes)
        --output_dir                       Path to output directory (must be surrounded with quotes)
        -profile                           Configuration profile to use. Can use multiple (comma separated)
    """.stripIndent()
}

params.help = false
params.hic_lib = false
params.graph = false

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

/*
if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check output_dir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.output_dir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}
*/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --              Miscellaneous code for the pipeline                    -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
import org.yaml.snakeyaml.Yaml

def parse_hic_lib_desc(filepath) {
    def yaml_file = new FileInputStream(new File(filepath))
    def lib_desc = new Yaml().load(yaml_file)
    params.hic_files = []
    lib_desc['fastqs'][0].eachWithIndex { entry, id ->
        params.hic_files.add([id, file(entry.get('R1'), checkIfExists: true), file(entry.get('R2'), checkIfExists: true)])
    }
    params.res_sites = lib_desc['res_sites'][0]
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       VALIDATE INPUTS                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if (params.help) {
    helpMessage()
    exit 0
}

if (!params.hic_lib || !params.graph) {
    exit 1, "Parameters '--hic_lib' and '--graph' are required to run the pipeline"
}

file(params.hic_lib, checkIfExists: true)
parse_hic_lib_desc(params.hic_lib)

file(params.graph, checkIfExists: true)

if (params.res_sites[0] == 'DNASE') {
    params.bwaopt = ''
} else {
    params.bwaopt = 'M'
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       HEADER LOG INFO                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// TODO

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                  SET UP INITIAL CHANNELS                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

ORIGINAL_GRAPH = Channel.fromPath(params.graph)
HIC_READS = Channel.fromList(params.hic_files)

// TODO

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          MAIN WORKFLOW                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

process preprocess_graph {
    publishDir "${params.output_dir}/", mode: 'copy'
    executor 'local'
    input:
        file ORIGINAL_GRAPH
    output:
        file 'fastas/bubbly_utg/' into fasta_files_for_bwa
        file 'fastas/sep_bubbly_utg/' into fasta_files_for_minimap2
        file 'fastas/non_bubbly_utg/'
        file 'inp_graph/'
        file 'logs/preprocess.log'
    script:
    """
        preprocess.py -g ${ORIGINAL_GRAPH} -o .
    """
}

process build_bwa_index {
    publishDir "${params.output_dir}/fastas/", mode: 'copy'
    executor 'local'
    input:
        file fasta_dir from fasta_files_for_bwa
    output:
        file fasta_dir into indexed_fasta_for_bwa
    script:
    """
         bwa index ${fasta_dir}/seqs.fa
    """
}

//mkdir ${fasta_file.getSimpleName()}_index && mv ${fasta_file}* ${fasta_file.getSimpleName()}_index

// --MD help for sniffles

process do_minimap2_alignment {
    publishDir "${params.output_dir}/pafs/", mode: 'copy'
    executor 'local'
    input:
        file fasta_dir from fasta_files_for_minimap2
    output:
        file "${fasta_dir[0].getSimpleName()}.paf.gz" into bubbly_alignment
    script:
    """
         minimap2 -xasm5 --cs=short ${fasta_dir}/seqs1.fa ${fasta_dir}/seqs2.fa | gzip -c - > ${fasta_dir[0].getSimpleName()}.paf.gz
    """
}


process do_bwa_alignment {
    maxForks 10
    publishDir "${params.output_dir}/bams/", mode: 'copy'
    echo true
    input:
        each file(index_dir) from indexed_fasta_for_bwa
        tuple val(id), file(reads1), file(reads2) from HIC_READS
    output:
        file "algn${id}.bam" into bamlets
    script:
    """
        bwa mem -SP5 -t 8 ${index_dir}/seqs.fa ${reads1} ${reads2} | samtools view -bS - > algn${id}.bam
    """
}



workflow.onComplete {
	log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

