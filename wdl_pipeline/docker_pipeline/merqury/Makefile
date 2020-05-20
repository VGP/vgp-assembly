# Definitions
repository = tpesout
identifier = vgp_merqury
version = 0.0.0
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build -t ${name}:${tag} ../../../dx_applets/merqury/docker/
	docker tag ${name}:${tag} ${name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest
