#!/bin/bash

corpus_folder=$1

echo ${corpus_folder}

now=$(date)
folder_name=results_$(date +'%s')
#cp -r template ${folder_name}
mkdir ${folder_name}

# Run our sequential block_tree construction algorithms
file_name=${folder_name}/repetitive_our.txt
lscpu | tee -a ${file_name}

for file in ${corpus_folder}*
do
    if [ -d ${file} ]; then
        continue
    fi
    for rep in 1 2 3
    do
    	for max_leaf_size in 2 4 8 16
    	do
            for tau in 2 4 8 16
            do
		for s in 0 1
		do
		    ./build/benchmark ${file} -e -c -s ${s} -t ${tau} -m ${max_leaf_size} -q 0 -p 1KB | tee -a ${file_name}
		done
	    done
        done
    done
done

# Run Belazzougui et al.'s sequential block_tree construction algorithm
file_name=${folder_name}/repetitive_belazzougui.txt
lscpu | tee -a ${file_name}

for file in ${corpus_folder}*
do
    if [ -d ${file} ]; then
        continue
    fi
    for rep in 1 2 3
    do
    	for max_leaf_size in 2 4 8 16
    	do
            for tau in 2 4 8 16
            do
		./build/benchmark ${file} -e -c -s 1 -t ${tau} -m ${max_leaf_size} -q 0 -p 1KB | tee -a ${file_name}
	    done
        done
    done
done

exit

# Run our parallel block tree construction algorithm
file_name=${folder_name}/repetitive_our_parallel.txt
lscpu | tee -a ${file_name}

for file in ${corpus_folder}*
do
    if [ -d ${file} ]; then
        continue
    fi
    for rep in 1 2 3
    do
    	for max_leaf_size in 2 4 8 16
    	do
            for tau in 2 4 8 16
            do
	        for threads in 1 2 4 8 16 32 64 128
	    	do
		    export OMP_NUM_THREADS=${threads}
		    ./build/benchmark ${file} -e -c -s 1 -t ${tau} -m ${max_leaf_size} -q 0 --threads ${threads} --parallel | tee -a ${file_name}
	    	done
	    done
        done
    done
done
