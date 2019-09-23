# Definitions
name = tpesout/modules
version = 4.3.0
git_commit = c999abdc958fe45c735dca093555b5413aecb678
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build -t ${name}:${tag} --build-arg git_commit=${git_commit} .
	docker tag ${name}:${tag} ${name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest
