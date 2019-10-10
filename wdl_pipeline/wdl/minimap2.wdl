workflow helloMinimap2 {
	call minimap2 
}

task minimap2 {
	File refFasta
	File readFile
	String dockerMinimap2
	command {
		cp ${refFasta} .
		cp ${readFile} .
		echo `basename ${readFile}` >>input.fofn
		echo `basename ${refFasta} | sed 's/.fasta$//g' | sed 's/.fa$//g'` >>outputBase
		docker run --rm -i -v `pwd`:/data ${dockerMinimap2} `basename ${refFasta}`
	}
	output {
		String outputBase = read_string("outputBase")
		File minimap2Bam = outputBase + ".bam"
		File minimap2BamIdx = outputBase + ".bam.bai"
	}
}
