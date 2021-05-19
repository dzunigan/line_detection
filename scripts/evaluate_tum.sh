#! /bin/bash

for MIN_LINE_PERCENT in 01 05 10 15 20 25 30 35 40 45 50
do
	./build/evaluate_tum --min_line=$MIN_LINE_PERCENT --output_file="./results_${MIN_LINE_PERCENT}.csv"
done
