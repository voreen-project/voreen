#!/usr/bin/env bash
set -e

if [ $# -ne 1 ]; then
    echo "usage: gen_compile_commands_json.sh <build_dir>"
    exit 1
fi

build_dir=$1

if command -v compdb > /dev/null ; then
    mkdir -p $build_dir

    cd $build_dir
    src_root="./$(git rev-parse --show-cdup)/voreen"
    cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 $src_root
    cd -

    output="./$(git rev-parse --show-cdup)/compile_commands.json"
    if compdb -p $build_dir list > $output; then
        echo "Successfully generated $output"
    else
        echo "Failed to generate $output"
    fi
else
    echo "compdb (https://github.com/Sarcasm/compdb) is not installed"
fi
