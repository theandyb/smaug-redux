# smaug-redux

## What I've done?

1. Downloaded and prepared reference data using `download_ref_data.sh` and `prep_ref_data.sh`

1. Annotated VCFs with motif and mutation type using `data_mgmt/data_prep/vcf_to_summary.pl`

1. Run `perl data_mgmt\data_prep\count_motifs_batch.pl` to generate motif counts

1. Run `perl data_mgmt\data_prep\extract_sing.pl` to extract singletons to `singletons/`

1. Run `R/analysis.R`
   * Does all the analyses we can do before fitting genomic features models

1. Run `perl data_mgmt/per-site_dp/process_glf_master.pl` to pull per-site read depth from GLF files

1. Run `perl data_mgmt/logit_scripts/build_data_worker.pl` to generate input data

1. Run `perl data_mgmt/logit_scripts/runmod_batch.pl` to run logit models
