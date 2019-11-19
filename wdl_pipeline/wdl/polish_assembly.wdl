import "tasks/minimap2.wdl" as minimap2
import "tasks/marginPolish.wdl" as marginPolish

workflow PolishAssembly {
    File ASSEMBLY_FILE
    Array[File] READ_FILES
    String SAMPLE_NAME
    File MARGIN_POLISH_PARAMS
    String MINIMAP_PRESET="map-pb"
    String SAMTOOLS_FILTER="-F 0x904"
    Int THREAD_COUNT
    Int MEMORY_GB=8
    String DOCKER_REPOSITORY="tpesout"
    String DOCKER_TAG="latest"

    # actual work
	call minimap2.minimap2 as align {
	    input:
            refFasta=ASSEMBLY_FILE,
            readFiles=READ_FILES,
            minimapPreset=MINIMAP_PRESET,
            samtoolsFilter=SAMTOOLS_FILTER,
            dockerRepository=DOCKER_REPOSITORY,
            dockerTag=DOCKER_TAG
	}
	call marginPolish.marginPolish as polish {
	    input:
            sampleName=SAMPLE_NAME,
            alignmentBam=align.minimap2Bam,
            alignmentBamIdx=align.minimap2BamIdx,
            referenceFasta=ASSEMBLY_FILE,
            parameters=MARGIN_POLISH_PARAMS,
            featureType="",
            threadCount=THREAD_COUNT,
            memoryGigabyte=MEMORY_GB,
            dockerRepository=DOCKER_REPOSITORY,
            dockerTag=DOCKER_TAG
	}

	output {
		File polishedAssembly = polish.polishedFasta
	}
}
