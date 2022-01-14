# 40K-RNA

Git repo for amalgamating splicing events from various public (and non-public) RNASeq datasets. Developed in R version 3.6 and human genome version hg19.

| Data source | URL                                   | Notes                                                        | Where to place the file?            |
| ----------- | ------------------------------------- | ------------------------------------------------------------ | ----------------------------------- |
| Intropolis  | https://github.com/nellore/intropolis | Based on  [intropolis.v1.hg19.tsv.gz](http://bit.ly/1SfBRTi)<br /> | Intropolis_sj_stats / data / source |
| GTEx v8     | dbGap V8                              | Individual SJ output downloaded for all GTEx samples<br /><br />These are hg38. Post processing scripts  combined these into a single tab file & lift over to hg19 | gtex_sj_stats / data / processed    |

# output

pre-processed 40K-RNA is available [here](https://storage.googleapis.com/misspl-db-data/misspl_events_40k_hg19.sql.gz)

# Instructions
1. `gtex_sj_stats` contains code to combine individual GTEx SJ files and summarise. 
first run:
  `snakemake -d . --cores 4 all`
then run:
  `python transform_bed_to_tsv.py data/processed/GTEx_V8.liftover_hg19.sorted.bed.gz   data/processed/GTEx_V8.liftover_hg19.sorted.tsv.gz`

2. `intropolis_sj_stats` contains code to summarise intropolis with max reads across samples for each SJ
Run:
  `Rscript calculate_intropolis_sj_summary.R "data/source/intropolis.v1.hg19.tsv.gz" "data/processed/intropolis_SJ_summary.sqlite" 20000`
  
3. R notebooks process GTEx and Intropolis data further, processing splice-junctions and infers which mis-splicing events they correspond to for each annotated donor & acceptor. The end result is 40K-RNA.

Run: `gtex_processing.Rmd` and `intropolis_processing.Rmd` 
Then run: `get_predictions.Rmd`. This will correlate observed SJ with mis-splicing events at each ensembl-annotated splice site.

