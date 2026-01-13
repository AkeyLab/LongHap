This Snakemake workflow reproduces the benchmarking of LongHap and other read-based variant phasing methods (WhatsHap, HapCUT2, LongPhase, and WhatsHap + MethPhaser) on HG002 for Pfennig and Akey (2025).

It assumes that the PacBio HiFi, ONT, and UL-ONT data are already downloaded. And stored in the specified `{data_dir}/HG002/` directory. The HiFi data (bam) should be stored in a subdirectory called `hifi_reads/`, the ONT data (bam) in a subdirectory called `ul_ont_r10.4.1_dorado/`, and the UL-ONT data (fastqs) in a subdirectory called `ul_ont_r10.4.1_dorado/`.
We retrieved the sequencing data from the following sources:
- PacBio HiFi: https://human-pangenomics.s3.amazonaws.com/submissions/80d00e88-7a92-46d8-88c7-48f1486e11ed--HG002_PACBIO_REVIO/m84039_230117_233243_s1.hifi_reads.default.bam
- ONT R10.4.1 Dorado base calling: s3://ont-open-data/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam
- UL-ONT R10.4.1 Dorado base calling: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/HG002/sequencing/ont/11_16_22_R1041_UL_HG002_dorado0.4.0/


The workflow will download the T2T reference genome (hs1) and the HG002 T2T-Q100 benchmark variant calls for HG002.

To run the workflow, simply execute:
```commandline
snakemake --profile profiles/ --use-conda
```

The figures can then be reproduced using Jupyter Notebook (`benchmark_longhap.ipynb`).