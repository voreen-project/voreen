#!/bin/bash
mkdir empty
docker build -f DockerfileBuildRunner -t "voreen/ci-ubuntu-16.04" empty
rmdir empty
