# smaug-redux

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
for i in \`seq 1 22\`; do
python data_mgmt/data_prep/motif_count.py -i /path/to/smaug-redux/reference_data/human_g1k_v37/chr$i.fasta -m data_mgmt/data_prep/motifs7.txt -o /path/to/smaug-redux/motif_counts/7-mers/full -c $i -b /path/to/smaug-redux/reference_data/genome.1000kb.sorted.bed
done
```

4. `perl data_mgmt/data_prep/extract_sing.pl` which runs vcftools --singletons on the VCF files, outputing results to `singletons/` directory

5. First pass through `analysis.r`

## To-do / Wishlist

1. Skip `extract_sing.pl` step -> get subject ids when generating summary files

2. Use `motif_count.py` for all motif counting

    * Generate counts in bins and then sum?
    
    * Second option for getting the counts across the entire genome?

3. Run the logistic regression analyses

4. Run the model evaluations 
