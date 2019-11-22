import "tasks/minimap2_scatter_gather.wdl" as minimap2
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
	call minimap2.runMinimap2ScatterGather as align {
	    input:
            refFasta=ASSEMBLY_FILE,
            readFiles=READ_FILES,
            minimapPreset=MINIMAP_PRESET,
            samtoolsFilter=SAMTOOLS_FILTER,
            threadCount=THREAD_COUNT,
            dockerRepository=DOCKER_REPOSITORY,
            dockerTag=DOCKER_TAG
	}
	call marginPolish.marginPolish as polish {
	    input:
            sampleName=SAMPLE_NAME,
            alignmentBam=align.alignment,
            alignmentBamIdx=align.alignmentIdx,
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


#import "tasks/minimap2.wdl" as minimap2
#import "tasks/marginPolish.wdl" as marginPolish
#
#workflow PolishAssembly {
#    File ASSEMBLY_FILE
#    Array[File] READ_FILES
#    String SAMPLE_NAME
#    File MARGIN_POLISH_PARAMS
#    String MINIMAP_PRESET="map-pb"
#    String SAMTOOLS_FILTER="-F 0x904"
#    Int THREAD_COUNT
#    Int MEMORY_GB=8
#    String DOCKER_REPOSITORY="tpesout"
#    String DOCKER_TAG="latest"
#
#    # actual work
#	call minimap2.minimap2_idx as idx {
#        input:
#            refFasta=ASSEMBLY_FILE,
#            threadCount=THREAD_COUNT,
#            minimapPreset=MINIMAP_PRESET,
#            dockerRepository=DOCKER_REPOSITORY,
#            dockerTag=DOCKER_TAG,
#	}
#    scatter (file in READ_FILES) {
#        call minimap2.minimap2_align as aln {
#            input:
#                refFasta=ASSEMBLY_FILE,
#                refFastaIdx=idx.minimap2Index,
#                readFile=file,
#                threadCount=THREAD_COUNT,
#                minimapPreset=MINIMAP_PRESET,
#                samtoolsFilter=SAMTOOLS_FILTER,
#                dockerRepository=DOCKER_REPOSITORY,
#                dockerTag=DOCKER_TAG
#            }
#    }
#	call minimap2.minimap2_merge as merge {
#        input:
#            outputBase=idx.outputBase,
#            alignFiles=aln.minimap2Bam,
#            alignIndexFiles=aln.minimap2BamIdx,
#            threadCount=THREAD_COUNT,
#            dockerRepository=DOCKER_REPOSITORY,
#            dockerTag=DOCKER_TAG,
#	}
#	call marginPolish.marginPolish as polish {
#	    input:
#            sampleName=SAMPLE_NAME,
#            alignmentBam=merge.alignment,
#            alignmentBamIdx=merge.alignmentIdx,
#            referenceFasta=ASSEMBLY_FILE,
#            parameters=MARGIN_POLISH_PARAMS,
#            featureType="",
#            threadCount=THREAD_COUNT,
#            memoryGigabyte=MEMORY_GB,
#            dockerRepository=DOCKER_REPOSITORY,
#            dockerTag=DOCKER_TAG
#	}
#	output {
#		File polishedAssembly = polish.polishedFasta
#	}
#}
