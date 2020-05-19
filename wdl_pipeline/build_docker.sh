#!/usr/bin/env bash
set -e

export REPOSITORY="tpesout"
if [[ ! -z "$1" ]] ; then export REPOSITORY="$1" ; fi

export VERSION="0.1.0"
if [[ ! -z "$2" ]] ; then export VERSION="$2" ; fi

export TAG_KEY="version"
if [[ "$3" == "release" ]] || [[ "$3" == "RELEASE" ]] ; then export TAG_KEY="tag" ; fi

echo "Building and deploying all VGP Docker images with:"
echo "  repository: '$REPOSITORY'"
echo "  version:    '$VERSION'"
if [[ $TAG_KEY == "tag" ]] ; then echo "  RELEASE" ; fi
sleep 2

# base
cd docker_base
make repository=$REPOSITORY $TAG_KEY=$VERSION
cd ..

# build scripts
cd docker_pipeline
for DIR in `ls | grep --invert-match bionano` ; do
    cd $DIR
    make repository=$REPOSITORY $TAG_KEY=$VERSION
    cd ..
done

# push all
cd ../docker_base
make push repository=$REPOSITORY version=$VERSION
cd ../docker_pipeline
for DIR in `ls | grep --invert-match bionano` ; do
    cd $DIR
    make push repository=$REPOSITORY $TAG_KEY=$VERSION
    cd ..
done

echo
echo "Remember, bionano still needs to be built manually."
echo
echo "Fin."