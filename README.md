# mtscATAC_pilot_analysis
##pilot data and analysis scripts for an mtscATAC-seq pilot performed on aged mice brain tissue.
##Daniel Richard 2024

The relevant data files (mostly the fragments.tsv.gz, cell metadata and mgatk output) are here: https://www.dropbox.com/scl/fo/mlx87ar4013f1tbkvj6wp/AID2i1Ziel7fQPcTip2FKnw?rlkey=67twjo25fzbuhh62dpz0hdnrz&st=l528t1uf&dl=0

Most analysis revolves around the ArchR script. There's a secondary script for how the 10X data was loaded into R for eventual use with FindClonotypes.

MGATK was run using the latest version available on their Github. All defaults were used, with the command

mgatk tenx  -g mm10 -i possorted_bam.bam -n $SAMPLE -o $OUTDIR -c 24 -bt CB -b barcodes.tsv --keep-temp-files --keep-qc-bams
