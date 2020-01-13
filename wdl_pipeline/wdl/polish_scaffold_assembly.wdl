version 1.0

import "tasks/extract_reads.wdl" as extractReads_t
import "tasks/shasta.wdl" as shasta_t
import "tasks/minimap2_scatter_gather.wdl" as minimap2_t
import "tasks/marginPolish.wdl" as marginPolish_t
import "tasks/purge_dups.wdl" as purgeDups_t
import "tasks/scaff10x.wdl" as scaff10x_t
import "tasks/salsa.wdl" as salsa_t
import "tasks/busco.wdl" as busco_t
import "tasks/stats.wdl" as stats_t

workflow ONTAssembly {
    input {
        File ASSEMBLY_FASTA
        Array[File] READ_FILES_POLISH
        Array[File] READ_FILES_10X
        Array[File] READ_FILES_HIC
        String SAMPLE_NAME
        File MARGIN_POLISH_PARAMS
        Int EXPECTED_GENOME_SIZE
        Int THREAD_COUNT
        Int? MEMORY_GB=8
        String? DOCKER_REPOSITORY="tpesout"
        String? DOCKER_TAG="latest"
    }

    Int default_MEMORY_GB = select_first([MEMORY_GB, 32])
    String default_DOCKER_REPOSITORY = select_first([DOCKER_REPOSITORY, "tpesout"])
    String default_DOCKER_TAG = select_first([DOCKER_TAG, "latest"])

    # actual work
    scatter (readFile in READ_FILES_POLISH) {
        call extractReads_t.extractReads as polishReads {
            input:
                readFile=readFile,
                dockerRepository=default_DOCKER_REPOSITORY,
                dockerTag=default_DOCKER_TAG
        }
    }

	call minimap2_t.runMinimap2ScatterGather as assemblyAlign {
	    input:
            refFasta=ASSEMBLY_FASTA,
            readFiles=polishReads.outputFile,
            minimapPreset="map-ont",
            samtoolsFilter="-F 0x904",
            threadCount=THREAD_COUNT,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call marginPolish_t.marginPolish as marginPolish {
	    input:
            sampleName=SAMPLE_NAME,
            alignmentBam=assemblyAlign.alignment,
            alignmentBamIdx=assemblyAlign.alignmentIdx,
            referenceFasta=ASSEMBLY_FASTA,
            parameters=MARGIN_POLISH_PARAMS,
            featureType="",
            threadCount=THREAD_COUNT,
            memoryGigabyte=default_MEMORY_GB,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
    call purgeDups_t.purge_dups as purgeDups {
        input:
            assemblyFasta=marginPolish.polishedFasta,
            readFiles=polishReads.outputFile,
            minimapPreset="map-ont",
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            memoryGigabyte=default_MEMORY_GB,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
    }
    call scaff10x_t.scaff10x as scaff10x {
        input:
            assemblyFasta=purgeDups.primary,
            readFiles10x=READ_FILES_10X,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            memoryGigabyte=default_MEMORY_GB,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
    }
    call salsa_t.salsa as salsa {
        input:
            refFasta=scaff10x.scaffoldedFasta,
            readFilesHiC=READ_FILES_HIC,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
    }
    

	### stats and validation
	# asm stats
	call stats_t.stats as assembly_stats {
	    input:
	        assemblyFasta=ASSEMBLY_FASTA,
	        expectedGenomeSize=EXPECTED_GENOME_SIZE,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call stats_t.stats as marginPolish_stats {
	    input:
	        assemblyFasta=marginPolish.polishedFasta,
	        expectedGenomeSize=EXPECTED_GENOME_SIZE,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call stats_t.stats as purgeDupsPri_stats {
	    input:
	        assemblyFasta=purgeDups.primary,
	        expectedGenomeSize=EXPECTED_GENOME_SIZE,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call stats_t.stats as purgeDupsAlt_stats {
	    input:
	        assemblyFasta=purgeDups.alternate,
	        expectedGenomeSize=EXPECTED_GENOME_SIZE,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call stats_t.stats as scaff10x_stats {
	    input:
	        assemblyFasta=scaff10x.scaffoldedFasta,
	        expectedGenomeSize=EXPECTED_GENOME_SIZE,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call stats_t.stats as salsa_stats {
	    input:
	        assemblyFasta=salsa.scaffoldedFasta,
	        expectedGenomeSize=EXPECTED_GENOME_SIZE,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	# busco stats
	call busco_t.busco as assembly_busco {
	    input:
	        assemblyFasta=ASSEMBLY_FASTA,
            threadCount=THREAD_COUNT,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call busco_t.busco as marginPolish_busco {
	    input:
	        assemblyFasta=marginPolish.polishedFasta,
            threadCount=THREAD_COUNT,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call busco_t.busco as purgeDupsPrimary_busco {
	    input:
	        assemblyFasta=purgeDups.primary,
            threadCount=THREAD_COUNT,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call busco_t.busco as scaff10x_busco {
	    input:
	        assemblyFasta=scaff10x.scaffoldedFasta,
            threadCount=THREAD_COUNT,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}
	call busco_t.busco as salsa_busco {
	    input:
	        assemblyFasta=salsa.scaffoldedFasta,
            threadCount=THREAD_COUNT,
            dockerRepository=default_DOCKER_REPOSITORY,
            dockerTag=default_DOCKER_TAG
	}

	output {
		File marginPolishFasta = marginPolish.polishedFasta
		File purgedPrimaryFasta = purgeDups.primary
		File purgedAlternateFasta = purgeDups.alternate
		File scaff10xFasta = scaff10x.scaffoldedFasta
		File salsaFasta = salsa.scaffoldedFasta

        # validation
		File assemblyBuscoResult = assembly_busco.outputTarball
		File marginPolishBuscoResult = marginPolish_busco.outputTarball
		File purgeDupsPrimaryBuscoResult = purgeDupsPrimary_busco.outputTarball
		File scaff10xBuscoResult = scaff10x_busco.outputTarball
		File salsaBuscoResult = salsa_busco.outputTarball

		File assemblyStatsResult = assembly_stats.statsTarball
		File marginPolishStatsResult = marginPolish_stats.statsTarball
		File purgedPrimaryStatsResult = purgeDupsPri_stats.statsTarball
		File purgedAlternateStatsResult = purgeDupsAlt_stats.statsTarball
		File scaff10xStatsResult = scaff10x_stats.statsTarball
		File salsaStatsResult = salsa_stats.statsTarball
	}
}
