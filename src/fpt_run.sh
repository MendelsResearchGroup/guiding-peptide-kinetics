#!/usr/bin/env bash

for i in {0..4}; do
    ./src/fpt_single_run.sh "$i" --force &
    [[ $(jobs -r -p | wc -l) -ge 5 ]] && wait -n
done
wait
