version 1.0

import "https://raw.githubusercontent.com/tpesout/vgp-assembly/b18a054e614d8bbfdf4d34936fc36b5174356a47/wdl_pipeline/wdl/tasks/extract_reads.wdl" as extractReads
import "https://raw.githubusercontent.com/tpesout/vgp-assembly/b18a054e614d8bbfdf4d34936fc36b5174356a47/wdl_pipeline/wdl/tasks/minimap2.wdl" as minimap2
import "https://raw.githubusercontent.com/tpesout/vgp-assembly/b18a054e614d8bbfdf4d34936fc36b5174356a47/wdl_pipeline/wdl/tasks/marginPolish.wdl" as marginPolish

workflow PolishAssembly {
    input {
        File ASSEMBLY_FILE
        Array[File] READ_FILES
        String SAMPLE_NAME
        File MARGIN_POLISH_PARAMS
        String? MINIMAP_PRESET="map-ont"
        String? SAMTOOLS_FILTER="-F 0x904"
        Int THREAD_COUNT
        Int? MARGINPOLISH_MEMORY_GB=32
        String? DOCKER_REPOSITORY="tpesout"
        String? DOCKER_TAG="latest"
    }

    call extractReads.runExtractReads as extract {
        input:
            inputFiles=READ_FILES,
            dockerImage=DOCKER_REPOSITORY+"/vgp_base:"+DOCKER_TAG
    }

    call minimap2.runMinimap2ScatterGather as align {
        input:
            refFasta=ASSEMBLY_FILE,
            readFiles=extract.reads,
            minimapPreset=MINIMAP_PRESET,
            samtoolsFilter=SAMTOOLS_FILTER,
            threadCount=THREAD_COUNT,
            dockerImage=DOCKER_REPOSITORY+"/vgp_minimap2:"+DOCKER_TAG
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
            memoryGigabyte=MARGINPOLISH_MEMORY_GB,
            dockerImage=DOCKER_REPOSITORY+"/vgp_marginpolish:"+DOCKER_TAG
    }

    output {
        File polishedAssembly = polish.polishedFasta
    }

}
