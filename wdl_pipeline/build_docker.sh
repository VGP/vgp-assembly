#!/usr/bin/env bash
set -e

export REPOSITORY="tpesout"
if [[ ! -z "$1" ]] ; then
    export REPOSITORY="$1"
fi

export VERSION="0.0.1"
if [[ ! -z "$2" ]] ; then
    export VERSION="$2"
fi

echo "Building and deploying all VGP Docker images with:"
echo "  repository: '$REPOSITORY'"
echo "  version:    '$VERSION'"
sleep 2

# base
cd docker_base
make
cd ..

# build scripts
cd docker_pipeline
for DIR in `ls` ; do
    cd $DIR
    make repository=$REPOSITORY version=$VERSION
    cd ..
done

# push all
cd ../docker_base
make push repository=$REPOSITORY version=$VERSION
cd ../docker_pipeline
for DIR in `ls` ; do
    cd $DIR
    make push repository=$REPOSITORY version=$VERSION
    cd ..
done

echo "Fin."