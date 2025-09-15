Feature counts output:
---------------------------
gene1: 10 reads total
gene2: 8 reads total
gene3: 2 reads total
---------------------------


#test command
snakemake --cores 10 --use-singularity --keep-going --config \
processing_dir="tests/test_counts/test_out/processing" \
results_dir="tests/test_counts/test_out/results" \
benchmark_dir="tests/test_counts/test_out/benchmarks" \
genome_index="tests/test_counts/genome/random_genome" \
gtf_file="tests/test_counts/genome/annotation.gtf" \
samples_csv="test_samples.csv" \
cleanup_processing='True'
