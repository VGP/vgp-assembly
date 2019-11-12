import "shasta.wdl" as shasta
import "minimap2.wdl" as minimap2
import "marginPolish.wdl" as marginPolish
import "busco.wdl" as busco
import "stats.wdl" as stats

workflow ONTAssembly {
    Array[File] READ_FILES
    String SAMPLE_NAME
    File MARGIN_POLISH_PARAMS
    Int EXPECTED_GENOME_SIZE
    Int THREAD_COUNT
    Int MEMORY_GB


    call shasta.shasta as shastaAssemble {
        input:
            readFilesONT=READ_FILES,
            sampleName=SAMPLE_NAME,
            threadCount=THREAD_COUNT,
            memoryGigabyte=MEMORY_GB
    }
	call busco.busco as shastaBusco {
	    input:
	        assemblyFasta=shastaAssemble.assemblyFasta
	}
	call stats.stats as shastaStats {
	    input:
	        assemblyFasta=shastaAssemble.assemblyFasta,
	        expectedGenomeSize=EXPECTED_GENOME_SIZE
	}
	call minimap2.minimap2 as shastaAlign {
	    input:
            refFasta=shastaAssemble.assemblyFasta,
            readFiles=READ_FILES,
            minimapPreset="map-ont",
            samtoolsFilter="-F 0x904"
	}
	call marginPolish.marginPolish as shastaMarginPolish {
	    input:
            sampleName=SAMPLE_NAME,
            alignmentBam=shastaAlign.minimap2Bam,
            alignmentBamIdx=shastaAlign.minimap2BamIdx,
            referenceFasta=shastaAssemble.assemblyFasta,
            parameters=MARGIN_POLISH_PARAMS,
            featureType="",
            threadCount=THREAD_COUNT,
            memoryGigabyte=MEMORY_GB
	}
	call stats.stats as marginPolishStats {
	    input:
	        assemblyFasta=shastaMarginPolish.polishedFasta,
	        expectedGenomeSize=EXPECTED_GENOME_SIZE
	}
	call busco.busco as marginPolishBusco {
	    input:
	        assemblyFasta=shastaMarginPolish.polishedFasta
	}

	output {
		File shastaAssembly = shastaAssemble.assemblyFasta
		File shastaBuscoResult = shastaBusco.outputTarball
		File shastaStatsResult = shastaStats.statsTarball
		File marginPolishAssembly = shastaMarginPolish.polishedFasta
		File marginPolishBuscoResult = marginPolishBusco.outputTarball
		File marginPolishStatsResult = marginPolishStats.statsTarball
	}
}
