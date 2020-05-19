version 1.0

import "tasks/extract_reads.wdl" as extractReads_t
import "tasks/shasta.wdl" as shasta_t
import "tasks/purge_dups.wdl" as purgeDups_t
import "tasks/scaff10x.wdl" as scaff10x_t
import "tasks/bionano.wdl" as bionano_t
import "tasks/salsa.wdl" as salsa_t
import "tasks/busco.wdl" as busco_t
import "tasks/stats.wdl" as stats_t

workflow ScaffoldAssembly {
    input {
        File ASSEMBLY_FASTA
        Array[File] READ_FILES
        Array[File] READ_FILES_10X
        Array[File] CMAPS_FILES_BIONANO
        Array[File] READ_FILES_HIC
        String SAMPLE_NAME
        Int THREAD_COUNT
        Int? MEMORY_GB=32
        String? DOCKER_REPOSITORY="tpesout"
        String? DOCKER_TAG="latest"
    }

    # actual work
    scatter (readFile in READ_FILES) {
        call extractReads_t.extractReads as polishReads {
            input:
                readFile=readFile,
                dockerImage=DOCKER_REPOSITORY+"/vgp_base:"+DOCKER_TAG
        }
    }
    call purgeDups_t.purge_dups as purgeDups {
        input:
            assemblyFasta=ASSEMBLY_FASTA,
            readFiles=polishReads.outputFile,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            memoryGigabyte=MEMORY_GB,
            dockerImage=DOCKER_REPOSITORY+"/vgp_purge_dups:"+DOCKER_TAG
    }
    call scaff10x_t.scaff10x as scaff10x {
        input:
            assemblyFasta=purgeDups.primary,
            readFiles10x=READ_FILES_10X,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            memoryGigabyte=MEMORY_GB,
            dockerImage=DOCKER_REPOSITORY+"/vgp_scaff10x:"+DOCKER_TAG
    }
    call bionano_t.bionano as bionano {
        input:
            assemblyFasta=scaff10x.scaffoldedFasta,
            bionanoFiles=CMAPS_FILES_BIONANO,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            memoryGigabyte=MEMORY_GB,
            dockerImage="gcr.io/nanopore_dev/vgp_bionano:"+DOCKER_TAG
    }
    call salsa_t.salsa as salsa {
        input:
            refFasta=bionano.scaffoldedNTrimmedAsm,
            readFilesHiC=READ_FILES_HIC,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            dockerImage=DOCKER_REPOSITORY+"/vgp_salsa:"+DOCKER_TAG
    }

	output {
		File purgedPrimaryFasta = purgeDups.primary
		File purgedAlternateFasta = purgeDups.alternate
		File scaff10xFasta = scaff10x.scaffoldedFasta
		File bionanoFasta = bionano.scaffoldedNTrimmedAsm
		File salsaFasta = salsa.scaffoldedFasta
	}
}
