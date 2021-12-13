#!/bin/bash

Fastmerge -T${3} ${2}/${2} $(awk -F. -v path=${2} '{printf path"/"$1" "}' $1)
