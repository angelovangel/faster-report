#!/usr/bin/env bash

# use workdir to be the parent of the fastq dir
# collect arguments and execute faster-report.R in the container

usage="$(basename "$0") [-p fastqpath] [-r regex] [-t type] [-d rundate] [-f flowcell] [-s saveraw] [-u subsample] [-o outfile] [-h]

Execute faster-report in docker container, using workdir to be parent of provided fastqpath. The html report will be
saved there.
Options:
    -h  show this help text
    -p  (required) path to folder with fastq reads
    - all other options that you want to pass to faster-report.R"

while getopts :hp: flag
do
   case "${flag}" in
      h) echo "$usage"; exit;;
      p) path=${OPTARG};;
      :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
     #\?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   esac
done

# workdir to be the parent of the fastq dir, this script is all about this

docker run \
    --mount type=bind,src="$HOME",target="$HOME" \
    -w $(dirname $path) \
    aangeloo/faster-report \
    -p $path \
    "$@"
echo -e "Report saved in $(dirname $path)"