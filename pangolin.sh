#!/bin/bash
source /home/everett/miniconda3/etc/profile.d/conda.sh
conda activate pangolin
pangolin -t 25 -o summaries/allGenomes_90_5.pangolin summaries/allGenomes_90_5.fasta
