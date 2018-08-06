# smaug-redux

## What I've done?

1. Downloaded and prepared reference data using `download_ref_data.sh` and `prep_ref_data.sh`

1. Annotated VCFs with motif and mutation type using `data_mgmt/data_prep/vcf_to_summary.pl`

1. Run `perl data_mgmt\data_prep\count_motif_batch.pl` to generate motif counts

1. Run `perl data_mgmt\data_prep\extract_sing.pl` to extract singletons to `singletons/`
