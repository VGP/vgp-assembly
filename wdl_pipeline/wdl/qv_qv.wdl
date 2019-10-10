workflow helloMinimap2 {
	call qv
}

task qv {
	File bam
	File bamIdx
	String dockerQvQv
	command {
		cp ${bam} aligned.bam
		cp ${bamIdx} aligned.bam.bai
		echo `basename ${bam} | sed 's/.bam//'` >>outputBase
		docker run --rm -i -v `pwd`:/data ${dockerQvQv} `cat outputBase`
	}
	output {
		String outputBase = read_string("outputBase")
		File minimap2Bam = outputBase + ".bam"
		File minimap2BamIdx = outputBase + ".bam.bai"
	}
}
