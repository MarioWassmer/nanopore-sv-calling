#!/usr/bin/env nextflow

/*
=======================================================================================================================
                                        Oxford Nanopore long reads SV calling
=======================================================================================================================
-----------------------------------------------------------------------------------------------------------------------
@Homepage
https://github.com/MarioWassmer/clinVac.git
-----------------------------------------------------------------------------------------------------------------------
@Documentation
https://github.com/MarioWassmer/clinVac
-----------------------------------------------------------------------------------------------------------------------
@Repository
https://github.com/MarioWassmer/clinVac.git
-----------------------------------------------------------------------------------------------------------------------
*/



def helpMessage() {
    clinVarLogo()
    log.info"""
    
    Usage: 

    nextflow run main.nf -profile conda --reads /path/to/reads.fastq.gz --genome /path/to/genome.fasta

    
    Requirements
                                    All tools required to run this pipeline are available in the environment
                                    file. To use it set -profile conda when executing the pipeline.  

    Mandatory arguments
        --reads                     Path to nanopore long reads (no spaces,otherwise surround with quotes)
        --genome                    Path to the reference genome fasta

    Optional:
        --help                      Shows this help message with usage information
        --bam                       Path to a coordinate sorted bam file created with minimap2 or ngmlr from
                                    nanopore reads. This skips the alignment and sort processes.
                                    Important: Sniffles and Svim require the MD tag in the bam file
        --outDir                    Store pipeline output in this directory.
        -w                          Working directory. [Default: current folder]
        --skipNanoStat              Skip NanoStat execution on the sequencing reads

    Alignments
        --aligner                   Long read aligner to use: 'minimap2' [Default, recommended] or 'ngmlr'
        --publishSam                Additionally publish alignment sam file to output folder (default: false)

    Variant calling
      Sniffles
        --skipSniffles              Skips SNIFFLES variant calling [false].
        --snifflesNoGenotype        Don't perform genotype estimation [true].
        --snifflesNoCluster         Don't mark SVs that occur together [true].

      Svim
        --skipSvim                  Skips SVIM variant calling [false].
        --svimNoInsSequences        Don't report list of insertions sequences from supporting reads in INFO tag [true].
        --svimNoSequenceAlleles     Don't use nucleotide sequences for alleles of deletions and inversions [true].

    """.stripIndent()
}




/*
=======================================================================================================================
                                        Methods go here
=======================================================================================================================
*/


// show usage information
if (params.help) {
    helpMessage()
    exit 0
}

// Check alignment mode
if (params.aligner != "ngmlr" && params.aligner != "minimap2") {
    helpMessage()
    exit 0
}

// Check mandatory params
if ( !params.genome || !params.reads) {
    helpMessage()
    exit 0
} else {

    // open input files
    Channel.fromPath("$params.genome").into{ genomeAlignment; genomeIndexing; genomeSvim }
    //genome_file = file(params.genome)
    Channel.fromPath("$params.reads").into{ readsStats; readsAlignment }
    //reads_file = file(params.reads)
    //bamFile = Channel.fromPath("${params.bam}").ifEmpty{Channel.empty()}

/*
    // If pam file is specified
    if(params.bam) {
        // Index it if index is missing
        if (!file("${file(params.bam)}.bai").exists()) {
            process indexInputBam {
                publishDir("${file(params.bam).getParent()}", mode: "copy")
                input:
                file bam from file(params.bam)
                output:
                file "$bam" into bamSniffles, bamSvim
                file "${bam}.bai" into bamIndexSniffles, bamIndexSvim

                """
                samtools index $bam
                """
            }
        } else {
            // Bam file is already indexed, create channels
            Channel.fromPath("${params.bam}").into{ bamSniffles; bamSvim }
            Channel.fromPath(file("${file(params.bam)}.bai")).into{ bamIndexSniffles; bamIndexSvim }
        }
    }
*/
    // get path of genome file for publishing its index if missing
    genomeFilePath = "${file(params.genome).getParent()}"

}

/*
    If the reference genome index doesn't exist
    this process is executed before starting the pipeline.
*/
if (!file("${file(params.genome)}.fai").exists()) {
    process indexGenome {
    publishDir ("$genomeFilePath", mode: "copy")
    input:
    file genome from genomeIndexing

    output:
    file "${genome}.fai" into genome_index

    """
    samtools faidx $genome
    """
    }
}


process statistics {
    publishDir ("$params.outDir", mode: "copy")
    input:
    file reads from readsStats
    output:
    file "NanoStat/nanostat.log" into nanostat_output

    when:
    !params.skipNanoStat

    """
    NanoStat --fastq ${reads} --outdir NanoStat --name nanostat.log
    """
}

process getMinSupport {
    input:
    file stats from nanostat_output

    when:
    !params.skipNanoStat
    
    script:
    """
    cat $stats  | grep "Total bases" \
                | sed 's/\\s\\s*/ /g' \
                | cut -d ' ' -f3 \
                | cut -d '.' -f 1 \
                | sed 's/,//g' \
                | awk '{print \$1/3200000000}' \
                | cut -d '.' -f 1
    """
}



/*
    Execute this process only if 
    no bam file is provided.
*/
process alignment {
    publishDir ("$params.outDir", mode: "copy", enabled: params.publishSam)
    input:
    file genome from genomeAlignment
    file reads from readsAlignment

    output:
    file "${reads.simpleName}_${params.aligner}_alignment.sam" into alignmentFile
        
    when:
    !params.bam
    
    script:
        if (params.aligner == "ngmlr")
            """
            ngmlr -t ${task.cpus} -x ont --bam-fix -r ${genome} -q ${reads} \
            -o ${reads.simpleName}_${params.aligner}_alignment.sam
            """
        else if (params.aligner == "minimap2")
            """
            minimap2 -t ${task.cpus} -ax map-ont -z 600,200 --MD -L\
            -o ${reads.simpleName}_${params.aligner}_alignment.sam ${genome} ${reads}
            """
}

/*
    Execute this process only after the 
    alignment process, i.e. when no bam 
    input is given.
*/
process createBam {
    publishDir ("$params.outDir", mode: "copy")
    input:
    file alignment from alignmentFile

    output:
    file "${alignment.simpleName}.bam" into createdBam

    when:
    !params.bam

    """
    samtools sort -O BAM -o ${alignment.simpleName}.bam ${alignment}
    """
}

//if(!params.bam) {
    process indexBam {
        publishDir ("${params.outDir}", mode: "copy")
        input:
        file bam from createdBam

        output:
        file "$bam" into bamSniffles, bamSvim
        file "${bam}.bai" into bamIndexSniffles, bamIndexSvim

        //when:
        //!params.bam
        
        """
        samtools index $bam
        """
    }
//}

process sniffles {
    publishDir ("$params.outDir/unfiltered_calls/sniffles", mode: "copy")
    input:
    file bamIdx from bamIndexSniffles
    file bam from bamSniffles
        
    output:
    file "sniffles_results.vcf" into snifflesResultIntersect, snifflesResultUnion, snifflesEval
    file "sniffles.log" into snifflesLog

    when:
    !params.skipSniffles

    script:
    def minSupSniffles = params.snifflesMinSupport ? "--min_support ${params.snifflesMinSupport}" : "--min_support 7"
    def clusterSniffles = !params.snifflesNoCluster ? "--cluster" : ""
    def genotypeSniffles = !params.snifflesNoGenotype ? "--genotype" : ""
         
    """
    sniffles ${minSupSniffles} ${genotypeSniffles} ${clusterSniffles} -m ${bam} -v sniffles_results.vcf > sniffles.log
    """
}

process svim {
    publishDir ("$params.outDir/unfiltered_calls", mode: "copy")
    input:
    file bamIdx from bamIndexSvim
    file bam from bamSvim
    file genome from genomeSvim
    
    output:
    file "svim/*" into svim_all
    file "svim/final_results.vcf" into svimPostProcess, svimResultIntersect, svimResultUnion
    
    when:
    !params.skipSvim

    script:
    def insSeqSvim = !params.svimNoInsSequences ? "--insertion_sequences" : ""
    def seqAllelesSvim = !params.svimNoSequenceAlleles ? "--sequence_alleles" : ""

    """
    svim alignment ${insSeqSvim} ${seqAllelesSvim} svim ${bam} ${genome}
    """
}

/*
process filterSvimVCF {
    publishDir ("$params.outDir/final", mode: "copy")
    
    input:
    file svimInput from svimPostProcess
    
    output:
    file "svim_qc_filtered.vcf" into svim_prepared.vcf

    """
    cat <(cat ${svimInput} | grep "^#") \
        <(cat ${svimInput} | grep -vE "^#" \
        | egrep "^[(0-9)(X)(Y)(MT)]" \
        | awk -v threshold=11 '$6>=threshold {print $0}') \
        > svim_qc_filtered.vcf
    """
}

process highConfidence {
    publishDir("$params.outDir", mode: "copy")
    
    input:
    file sniffles from snifflesResultIntersect
    file svim from svimResultIntersect
    
    output:
    file "intersection.vcf" into highConfidenceSet
    
    // use survivor to get intersectio (get high confidence) of svim and sniffles
    """    
    SURVIVOR merge $IN 250 2 -1 -1 -1 49 intersection-set/intersection.vcf
    """
}

process highSensitivity {
    publishDir("$params.outDir", mode: "copy")

    input:
    file sniffles from snifflesResultUnion
    file svim from svimResultUnion

    output:
    file "intersection.vcf" into highSensitivitySet
    
    // use survivor to get union of svim and sniffles
    """
    """
}
*/


def clinVarLogo() {
    log.info " .----------------.  .----------------.  .----------------.  .-----------------. .----------------.  .----------------.  .----------------. ";
    log.info "| .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |";
    log.info "| |     ______   | || |   _____      | || |     _____    | || | ____  _____  | || | ____   ____  | || |      __      | || |     ______   | |";
    log.info "| |   .' ___  |  | || |  |_   _|     | || |    |_   _|   | || ||_   \\|_   _| | || ||_  _| |_  _| | || |     /  \\     | || |   .' ___  |  | |";
    log.info "| |  / .'   \\_|  | || |    | |       | || |      | |     | || |  |   \\ | |   | || |  \\ \\   / /   | || |    / /\\ \\    | || |  / .'   \\_|  | |";
    log.info "| |  | |         | || |    | |   _   | || |      | |     | || |  | |\\ \\| |   | || |   \\ \\ / /    | || |   / ____ \\   | || |  | |         | |";
    log.info "| |  \\ \\`.___.'\\ | || |   _| |__/ |  | || |     _| |_    | || | _| |_\\   |_  | || |    \\ ' /     | || | _/ /    \\ \\_ | || |  \\ \\`.___.'\\ | |";
    log.info "| |   \\`._____.' | || |  |________|  | || |    |_____|   | || ||_____|\\____| | || |     \\_/      | || ||____|  |____|| || |   \\`._____.' | |";
    log.info "| |              | || |              | || |              | || |              | || |              | || |              | || |              | |";
    log.info "| '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |";
    log.info " '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------' ";
}
