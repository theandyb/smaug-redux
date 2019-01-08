# smaug-redux

## Requirements

My suggestion is to run this analysis in a conda virtualenv

### R packages ()

r-tidyverse
r-devtools
r-broom 
r-rcolorbrewer
r-mass
r-speedglm
r-boot
r-devtools
r-psych
r-lmtest
r-fmsb
r-cowplot
r-hexbin
r-gridextra
r-gtable
r-nmf
bioconductor-biostrings
r-yaml
r-openxlsx
r-svglite
r-emdbook
r-bedr

### samtools-hybrid
Can be installed using conda with

```{bash}
conda install -c andyb262 samtools-hybrid 
```

## Instructions

1. Run `download_ref_data.sh` and `prep_ref_data.sh` to download and process necessary reference data

2. Run `perl data_mgmt/data_prep/vcf_to_summary.pl`
       
       * Creates a copy of the input vcf file, with additional annotations (ancestral allele, k-mer centered at that site, the mutation type)
        
        * Also generates summary files

3. `perl data_mgmt/data_prep/count_motifs_batch.pl`
        
        * Counts occurence of each motif in the reference genome
        
        * Use `data_mgmt/lib/motif_count.py` and `data_mgmt/lib/rc_motif_count.r` to get these counts in 1Mb windows

3. (option 2)

```bash
for i in `seq 1 22`; do
python data_mgmt/data_prep/motif_count.py -i /path/to/smaug-redux/reference_data/human_g1k_v37/chr$i.fasta -m data_mgmt/data_prep/motifs7.txt -o /path/to/smaug-redux/motif_counts/7-mers/full -c $i -b /path/to/smaug-redux/reference_data/genome.1000kb.sorted.bed
done
```

4. `perl data_mgmt/data_prep/extract_sing.pl` which runs vcftools --singletons on the VCF files, outputing results to `singletons/` directory

5. First pass through `analysis.r`

6. Follow the steps in `data_mgmt/per-site_dp/README.md` to extract depth information from glf files.

7. Run `perl data_mgmt/logit_scripts/build_data_worker.pl` to build input for regression models.

## To-do / Wishlist

1. Skip `extract_sing.pl` step -> get subject ids when generating summary files

2. Use `motif_count.py` for all motif counting

    * Generate counts in bins and then sum?
    
    * Second option for getting the counts across the entire genome?

3. Run the logistic regression analyses

4. Run the model evaluations 
