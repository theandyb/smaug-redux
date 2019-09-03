#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=1500MB
#SBATCH --ntasks=1
#SBATCH --time 00:30:00
#SBATCH --job-name=countMotifs
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue

python /net/snowwhite/home/beckandy/research/smaug-redux/data_mgmt/data_prep/motif_count.py -i /net/snowwhite/home/beckandy/research/smaug-redux/reference_data/human_g1k_v37/chr${SLURM_ARRAY_TASK_ID}.fasta -m /net/snowwhite/home/beckandy/research/smaug-redux/data_mgmt/data_prep/motifs7.txt -o /net/snowwhite/home/beckandy/research/smaug-redux/motif_counts/7-mers/full -c ${SLURM_ARRAY_TASK_ID}
