Feature counts output:
---------------------------
gene1: 10 reads total
gene2: 8 reads total
gene3: 2 reads total
---------------------------


snakemake --cores 10 --use-singularity --configfile test_config.yaml 

snakemake --cores 10 --use-singularity --config \
processing_dir="tests/test_counts/new_branch/processing" \
results_dir="tests/test_counts/new_branch/results" \
benchmark_dir="tests/test_counts/new_branch/benchmarks" \
genome_index="tests/test_counts/genome/random_genome" \
gtf_file="tests/test_counts/genome/annotation.gtf" \
samples_csv="test_samples_counts.csv"