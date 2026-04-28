### Benchmarking on HG002
This Snakemake workflow reproduces the benchmarking of LongHap and other read-based variant phasing methods (WhatsHap, HapCUT2, LongPhase, and WhatsHap + MethPhaser) on HG002 for Pfennig and Akey (2025).

It assumes that the PacBio HiFi, ONT, and UL-ONT data are already downloaded. And stored in the specified `{data_dir}/HG002/` directory. The HiFi data (bam) should be stored in a subdirectory called `hifi_reads/`, the ONT data (bam) in a subdirectory called `ul_ont_r10.4.1_dorado/`, and the UL-ONT data (fastqs) in a subdirectory called `ul_ont_r10.4.1_dorado/`.
We retrieved the sequencing data from the following sources:
- PacBio HiFi: https://human-pangenomics.s3.amazonaws.com/submissions/80d00e88-7a92-46d8-88c7-48f1486e11ed--HG002_PACBIO_REVIO/m84039_230117_233243_s1.hifi_reads.default.bam
- ONT R10.4.1 Dorado base calling: s3://ont-open-data/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam
- UL-ONT R10.4.1 Dorado base calling: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/HG002/sequencing/ont/11_16_22_R1041_UL_HG002_dorado0.4.0/


The workflow will download the T2T reference genome (hs1) and the HG002 T2T-Q100 benchmark variant calls for HG002.

You will have to download pb-CpG-tools and provide the path to the `aligned_bam_to_cpg_scores` executable in the `config/config.yaml`

The workflow also assumes that any aligner you're chosing to use (`minimap2`, `winnowmap2`, or `pbmm2`) is installed and available in your `PATH`.

To run the workflow, simply execute:
```commandline
snakemake -s Snakefile --profile profiles/ --use-conda
```

The figures can then be reproduced using Jupyter Notebook (`benchmark_longhap.ipynb` and `plot_locus.ipynb`).

### Population-level analysis
Additionally, we applied LongHap to 6 additional population-level samples using HiFi data from the HPRC:
- HG00609, CHS, (s3://human-pangenomics/working/HPRC/HG00609/raw_data/PacBio_HiFi/m84081_231124_014115_s1.hifi_reads.bc2090.bam, s3://human-pangenomics/working/HPRC/HG00609/raw_data/PacBio_HiFi/m84081_231124_021221_s2.hifi_reads.bc2090.bam, s3://human-pangenomics/working/HPRC/HG00609/raw_data/PacBio_HiFi/m84081_231124_024327_s3.hifi_reads.bc2090.bam)
- HG00658, CHS (s3://human-pangenomics/working/HPRC/HG00658/raw_data/PacBio_HiFi/m84091_230905_192642_s1.hifi_reads.bc1018.bam, s3://human-pangenomics/working/HPRC/HG00658/raw_data/PacBio_HiFi/m84091_230905_195701_s2.hifi_reads.bc1018.bam, s3://human-pangenomics/working/HPRC/HG00658/raw_data/PacBio_HiFi/m84091_230905_202807_s3.hifi_reads.bc1018.bam)
- HG00738, PUR (s3://human-pangenomics/working/HPRC/HG00738/raw_data/PacBio_HiFi/m84081_231112_034048_s4.hifi_reads.bc2079.bam)
- HG01099, PUR (s3://human-pangenomics/working/HPRC/HG01099/raw_data/PacBio_HiFi/m84081_231124_014115_s1.hifi_reads.bc2092.bam, s3://human-pangenomics/working/HPRC/HG01099/raw_data/PacBio_HiFi/m84081_231124_021221_s2.hifi_reads.bc2092.bam, s3://human-pangenomics/working/HPRC/HG01099/raw_data/PacBio_HiFi/m84081_231124_024327_s3.hifi_reads.bc2092.bam)
- HG02723, GWD (s3://human-pangenomics/submissions/548dd68a-3e67-44a6-8b39-954b7a8eb835--HPRC_REVIO_EA_2023/bc2012-HG02723/m84036_230317_175945_s2.hifi_reads.bc2012.bam, s3://human-pangenomics/submissions/548dd68a-3e67-44a6-8b39-954b7a8eb835--HPRC_REVIO_EA_2023/bc2012-HG02723/m84039_230303_012244_s3.hifi_reads.bc2012.bam, s3://human-pangenomics/submissions/548dd68a-3e67-44a6-8b39-954b7a8eb835--HPRC_REVIO_EA_2023/bc2012-HG02723/m84039_230314_213047_s2.hifi_reads.bc2012.bam, s3://human-pangenomics/submissions/548dd68a-3e67-44a6-8b39-954b7a8eb835--HPRC_REVIO_EA_2023/bc2012-HG02723/m84039_230316_193003_s2.hifi_reads.bc2012.bam)
- HG02615, GWD (s3://human-pangenomics/working/HPRC/HG02615/raw_data/PacBio_HiFi/m84081_231105_031800_s4.hifi_reads.default.bam)

We retrieved the corresponding assemblies:
- HG002 (s3://human-pangenomics/working/HPRC_PLUS/HG002/assemblies/Q100/pansn/hg002v1.1.mat_MT.PanSN.fa.gz, s3://human-pangenomics/working/HPRC_PLUS/HG002/assemblies/Q100/pansn/hg002v1.1.pat.PanSN.fa.gz)
- HG00609 (s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00609/assemblies/freeze_2/HG00609_mat_hprc_r2_v1.0.1.fa.gz, s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00609/assemblies/freeze_2/HG00609_pat_hprc_r2_v1.0.1.fa.gz)
- HG00658 (s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00658/assemblies/freeze_2/HG00658_mat_hprc_r2_v1.0.1.fa.gz, s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00658/assemblies/freeze_2/HG00658_pat_hprc_r2_v1.0.1.fa.gz)
- HG00738 (s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00738/assemblies/freeze_2/HG00738_mat_hprc_r2_v1.0.1.fa.gz, s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00738/assemblies/freeze_2/HG00738_pat_hprc_r2_v1.0.1.fa.gz)
- HG01099 (s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG01099/assemblies/freeze_2/HG01099_mat_hprc_r2_v1.0.1.fa.gz, s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG01099/assemblies/freeze_2/HG01099_pat_hprc_r2_v1.0.1.fa.gz)
- HG02615 (s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG02615/assemblies/freeze_2/HG02615_mat_hprc_r2_v1.0.1.fa.gz, s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG02615/assemblies/freeze_2/HG02615_pat_hprc_r2_v1.0.1.fa.gz)
- HG02723 (s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG02723/assemblies/freeze_2/HG02723_mat_hprc_r2_v1.0.1.fa.gz, s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG02723/assemblies/freeze_2/HG02723_pat_hprc_r2_v1.0.1.fa.gz)

To analyze these samples, you have manually download the data, update the paths at the top of each Snakemake file, and then run:
```commandline
snakemake -s full_longhap_pipeline.smk --profile profiles/ --use-conda
```

and to analyze LongHap's phasing, run:
```commandline
snakemake -s analyze_longhap_phasing.smk --profile profiles/ --use-conda
```

The figures can be reproduced using Jupyter Notebook (`population_analyses.ipynb`).