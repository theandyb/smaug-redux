# smaug-redux

## Instructions

1. Run `download_ref_data.sh` and `prep_ref_data.sh` to download and process necessary reference data
2. `perl data_mgmt/data_prep/vcf_to_summary.pl'
        * Creates a copy of the input vcf file, with additional annotations (ancestral allele, k-mer centered at that site, the mutation type)
        * Also generates summary files
3. `perl data_mgmt/data_prep/count_motifs_batch.pl`
        * Counts occurence of each motif in the reference genome
        * Use `data_mgmt/lib/motif_count.py` and `data_mgmt/lib/rc_motif_count.r` to get these counts in 1Mb windows
4. `perl data_mgmt/data_prep/extract_sing.pl` which runs vcftools --singletons on the VCF files, outputing results to `singletons/` directory
5. First pass through `analysis.r`
