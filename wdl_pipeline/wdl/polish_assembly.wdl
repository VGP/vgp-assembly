version 1.0

import "tasks/extract_reads.wdl" as extractReads
import "tasks/minimap2_scatter_gather.wdl" as minimap2
import "tasks/marginPolish.wdl" as marginPolish

workflow PolishAssembly {
    input {
        File ASSEMBLY_FILE
        Array[File] READ_FILES
        String SAMPLE_NAME
        File MARGIN_POLISH_PARAMS
        String? MINIMAP_PRESET
        String? SAMTOOLS_FILTER
        Int THREAD_COUNT
        Int? MARGINPOLISH_MEMORY_GB
        String? DOCKER_REPOSITORY
        String? DOCKER_TAG
    }

    String defaultMinimapPreset = select_first([MINIMAP_PRESET, "map-ont"])
    String defaultSamtoolsFilter = select_first([SAMTOOLS_FILTER, "-F 0x904"])
    Int defaultMarginPolishMemoryGB = select_first([MARGINPOLISH_MEMORY_GB, 8])
    String defaultDockerRepository = select_first([DOCKER_REPOSITORY, "tpesout"])
    String defaultDockerTag = select_first([DOCKER_TAG, "latest"])

    call extractReads.runExtractReads as extract {
        input:
            inputFiles=READ_FILES,
            dockerRepository=defaultDockerRepository,
            dockerTag=defaultDockerTag
    }

    call minimap2.runMinimap2ScatterGather as align {
        input:
            refFasta=ASSEMBLY_FILE,
            readFiles=extract.reads,
            minimapPreset=defaultMinimapPreset,
            samtoolsFilter=defaultSamtoolsFilter,
            threadCount=THREAD_COUNT,
            dockerRepository=defaultDockerRepository,
            dockerTag=defaultDockerTag
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
            memoryGigabyte=defaultMarginPolishMemoryGB,
            dockerRepository=defaultDockerRepository,
            dockerTag=defaultDockerTag
    }

    output {
        File polishedAssembly = polish.polishedFasta
    }

}
