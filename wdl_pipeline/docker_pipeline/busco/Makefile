# Definitions
repository = tpesout
identifier = vgp_marginpolish
version = 0.0.0
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Definitions
name = tpesout/vgp_busco
version = 0.0.0
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	mkdir tmp
	cp ../../../pipeline/busco/* tmp
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest
	rm -r tmp

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest
