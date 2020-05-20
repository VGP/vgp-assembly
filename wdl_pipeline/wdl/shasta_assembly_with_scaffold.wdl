version 1.0

import "tasks/extract_reads.wdl" as extractReads_t
import "tasks/shasta.wdl" as shasta_t
import "tasks/minimap2.wdl" as minimap2_t
import "tasks/marginPolish.wdl" as marginPolish_t
import "tasks/purge_dups.wdl" as purgeDups_t
import "tasks/scaff10x.wdl" as scaff10x_t
import "tasks/salsa.wdl" as salsa_t
import "tasks/busco.wdl" as busco_t
import "tasks/stats.wdl" as stats_t

workflow ShastaAssembly {
    input {
        Array[File] READ_FILES_ONT
        Array[File] READ_FILES_10X
        Array[File] READ_FILES_HIC
        Array[File] CMAPS_FILES_BIONANO
        String SAMPLE_NAME
        File MARGIN_POLISH_PARAMS
        File SHASTA_PARAMS
        Float EXPECTED_GENOME_SIZE
        Int THREAD_COUNT
        Int MEMORY_GB=32
        String DOCKER_REPOSITORY="tpesout"
        String DOCKER_TAG="latest"
    }

    # actual work
    scatter (readFile in READ_FILES_ONT) {
        call extractReads_t.extractReads as ontReads {
            input:
                readFile=readFile,
                dockerImage=DOCKER_REPOSITORY+"/vgp_base:"+DOCKER_TAG
        }
    }

    call shasta_t.shasta as shastaAssemble {
        input:
            readFilesONT=ontReads.outputFile,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            shastaParameters=SHASTA_PARAMS,
            memoryGigabyte=MEMORY_GB,
            dockerImage=DOCKER_REPOSITORY+"/vgp_shasta:"+DOCKER_TAG
    }
	call minimap2_t.runMinimap2ScatterGather as assemblyAlign {
	    input:
            refFasta=shastaAssemble.assemblyFasta,
            readFiles=ontReads.outputFile,
            minimapPreset="map-ont",
            samtoolsFilter="-F 0x904",
            threadCount=THREAD_COUNT,
            dockerImage=DOCKER_REPOSITORY+"/vgp_minimap2:"+DOCKER_TAG
	}
	call marginPolish_t.marginPolish as marginPolish {
	    input:
            sampleName=SAMPLE_NAME,
            alignmentBam=assemblyAlign.alignment,
            alignmentBamIdx=assemblyAlign.alignmentIdx,
            referenceFasta=shastaAssemble.assemblyFasta,
            parameters=MARGIN_POLISH_PARAMS,
            featureType="",
            threadCount=THREAD_COUNT,
            memoryGigabyte=MEMORY_GB,
            dockerImage=DOCKER_REPOSITORY+"/vgp_marginpolish:"+DOCKER_TAG
	}
    call purgeDups_t.purge_dups as purgeDups {
        input:
            assemblyFasta=marginPolish.polishedFasta,
            readFiles=ontReads.outputFile,
            minimapPreset="map-ont",
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
            refFasta=scaff10x.scaffoldedFasta,
            bionanoFiles=CMAPS_FILES_BIONANO,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            memoryGigabyte=MEMORY_GB,
            dockerImage=DOCKER_REPOSITORY+"/vgp_bionano:"+DOCKER_TAG
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
		File shastaFasta = shastaAssemble.assemblyFasta
		File marginPolishFasta = marginPolish.polishedFasta
		File purgedHaplotypeFasta = purgeDups.primary
		File purgedAlternateFasta = purgeDups.alternate
		File scaff10xFasta = scaff10x.scaffoldedFasta
		File bionanoFasta = bionano.scaffoldedNTrimmedAsm
		File salsaFasta = salsa.scaffoldedFasta
	}
}
