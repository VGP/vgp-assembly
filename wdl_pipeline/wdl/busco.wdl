workflow helloBusco {
	call busco
}

task busco {
	File assemblyFasta
	String dockerBusco
	command {
		cp ${assemblyFasta} .
		basename ${assemblyFasta} | sed 's/.fasta$//g' | sed 's/.fa$//g' >outputBase
		docker run --rm -i -v `pwd`:/data ${dockerBusco} `basename ${assemblyFasta}`
	}
	output {
		String outputBase = read_string("outputBase")
		File outputFolder = outputBase
	}
}
