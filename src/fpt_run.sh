#!/usr/bin/env bash

for i in {0..49}; do
    ./src/fpt_single_run.sh "$i" YYDPETGTWE --force &
    [[ $(jobs -r -p | wc -l) -ge 5 ]] && wait -n
done
wait
