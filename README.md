# mtscATAC_pilot_analysis
##pilot data and analysis scripts for an mtscATAC-seq pilot performed on aged mice brain tissue.
##Daniel Richard 2024

The relevant data files (mostly the fragments.tsv.gz, cell metadata and mgatk output) are here: https://www.dropbox.com/scl/fo/mlx87ar4013f1tbkvj6wp/AID2i1Ziel7fQPcTip2FKnw?rlkey=67twjo25fzbuhh62dpz0hdnrz&st=l528t1uf&dl=0

Most analysis revolves around the ArchR script. There's a secondary script for how the 10X data was loaded into R for eventual use with FindClonotypes.

MGATK was run using the latest version available on their Github. All defaults were used, with the command

mgatk tenx  -g mm10 -i possorted_bam.bam -n $SAMPLE -o $OUTDIR -c 24 -bt CB -b barcodes.tsv --keep-temp-files --keep-qc-bams

Here's a visual of the strand correlation - vs - VMR plot output from MGATK:
![image](https://github.com/user-attachments/assets/10a7df24-4c3d-4cc6-9ba9-1e9b6f53add9)

Here are the clusters I defined using ArchR, as well as their predicted identities based on a reference dataset.
![Screenshot 2024-08-02 at 5 27 20 PM](https://github.com/user-attachments/assets/8fa3a8f0-0446-438f-a7d2-b1885752044b)

And after analyses, here's an example of a clonotype heatmap generated from this data and analysis workflow:

 ![image](https://github.com/user-attachments/assets/fd6eb1fe-5301-4004-8f1f-dabc73f87c33)

