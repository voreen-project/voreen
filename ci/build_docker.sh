#!/bin/bash
script=$0
script_dir=$(dirname $script)

mkdir empty
docker build -f $script_dir/DockerfileBuildRunner -t "voreen/ci-ubuntu-18.04" empty
rmdir empty
