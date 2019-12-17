# Definitions
repository = tpesout
identifier = vgp_qv
version = 0.0.0
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	mkdir tmp
	cp ../../../pipeline/qv/* tmp
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest
	rm -r tmp

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest
