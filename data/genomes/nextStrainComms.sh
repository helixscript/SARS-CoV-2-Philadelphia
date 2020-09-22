augur align \
  --nthreads 30 \
  --sequences data/sequences.fasta \
  --reference-sequence config/SARS-CoV-2_outgroup.gb \
  --output results/aligned.fasta \
  --fill-gaps;

augur tree \
  --nthreads 30 \
  --alignment results/aligned.fasta \
  --output results/tree_raw.nwk;

augur refine \
  --tree results/tree_raw.nwk \
  --alignment results/aligned.fasta \
  --metadata data/metadata.tsv \
  --output-tree results/tree.nwk \
  --output-node-data results/branch_lengths.json \
  --timetree \
  --coalescent opt \
  --date-confidence \
  --date-inference marginal \
  --clock-filter-iqd 4;

