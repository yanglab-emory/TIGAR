#!/usr/bin/bash

out=$1
block_path_pre=$2
n_blocks=$3

start_row=1
for i in `seq 0 $n_blocks`; do
	block_path=${block_path_pre}_${i}.txt
	
	# number the lines then add bgzipped lines to final output
	nl -nln -s$'\t' -v${start_row} ${block_path} \
	>> ${out}

	# get starting row number for next block
	start_row=$(( $(grep -c $ "${block_path}") + start_row ))

	# remove block input file
	rm -f ${block_path}
done

echo 'Bgzip/tabix-ing final output LD file.'
bgzip -f ${out} && \
tabix -f -s3 -b4 -e4 -S1 ${out}.gz
