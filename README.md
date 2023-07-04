# Benchmark Scripts for pasta::block_tree

[![DOI](https://zenodo.org/badge/661411780.svg)](https://zenodo.org/badge/latestdoi/661411780)

This repository contains scripts that help to reproduce the results of the paper "Faster Block Tree Construction".

To clone this repository and build the benchmark please use the following commands.
```
git clone https://github.com/pasta-toolbox/block_tree_experiments --recursive
cd block_tree_experiments
cmake -DCMAKE_BUILD_TYPE=Release -B build
cmake --build build -j
```
This will compile all code needed for the benchmark in the `build` folder, which is used in the `run_benchmarks.sh` script.
Next, we can run the script.

```
./run_benchmarks.sh <path to Pizza&Chili repetitive corpus>
```

This script will copy the [`template`](template/) folder and rename it to `results_CURRENTDATE`, where `CURRENTDATE` can be decoded using the `date` command (`date --date @CURRENTDATE`).
Within the folder, you find a LaTex-Documents that contains all layout and styling information needed to recreate the plots.
To make the plots based on the newly created data, we rely on [sqlplot-tools][].
This tool reads the newly created results files and passes all values directly into the LaTeX files.
Please refer to [sqlplot-tools]' website for information on how to install it.
If you have installed [sqlplot-tools] simply run `sqlplot-tools <file_name.tex>` and compile `<file_name.tex>` to see the final results.

*Warning:* if you are using a different text corpus for the experiments, you have to adapt the plots.
Otherwise, no data will be shown, as the name of the inputs are used to populate the plots.

[sqlplot-tools]: https://github.com/bingmann/sqlplot-tools/
