#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd ~/programming/scivis-contest-2022/voreen/build/bin/
./voreentool -platform minimal --useCaching false -w $DIR/preprocessing_similarity.vws --script $DIR/preprocessing_similarity.py
