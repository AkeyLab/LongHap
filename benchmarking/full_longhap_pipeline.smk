import glob

configfile: 'config/config.yaml'
longhap_path = '~/LongHap/longhap.py'
workdir: config['base_results_dir']
base_results_dir = config["base_results_dir"]
data_dir = config['data_dir']
ref_dir = config["ref_dir"]
aligner = config['aligner']
approx_sequenced_coverage = config['approx_sequenced_coverage']
haplomaker_jar = config['haplomaker_jar']
aligned_bam_to_cpg_scores_path = config['aligned_bam_to_cpg_scores_path']
reference_genome_url = config['reference_genome_url']
par_regions_bed_url = config['par_regions_bed_url']
gt_variants_url = config['gt_variants_url']
autosomes = [str(i) for i in range(1, 23)]
sex_chromosomes=['X', "Y"]
phasing_methods = ['longhap_meth', "longhap"]

samples = ['HG00609', 'HG00658', 'HG00738', 'HG01099', 'HG02723', 'HG02615', 'HG002']

seq_tech = ['hifi']

rule all:
    input:
        expand(base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.vcf.gz.tbi",
               chrom=autosomes,  phasing_method=phasing_methods, aligner=aligner, seq_tech=seq_tech, sample=samples),
        expand(base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.stats.tsv",
               chrom=autosomes, phasing_method=phasing_methods, aligner=aligner, seq_tech=seq_tech, sample=samples),
        expand(base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.stats_complex_variants.tsv",
               chrom=autosomes, phasing_method=phasing_methods, aligner=aligner, seq_tech=seq_tech, sample=samples),

rule download_reference_genome:
    output:
        ref_dir + reference_genome_url.split('/')[-1]
    params:
        url = reference_genome_url,
        outdir = ref_dir
    conda:
        "envs/samtools.yaml"
    localrule: True
    shell:
        "mkdir -p {params.outdir}; wget -qO - {params.url} | gunzip | bgzip > {output}"
    
rule unzip_reference_fasta:
    input:
        ref=ref_dir + reference_genome_url.split('/')[-1]
    output:
        unzip=ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".fa")
    localrule: True
    shell:
        "gunzip -c {input.ref} > {output.unzip}; "
    
rule reference_to_upper:
    input:
        ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', '.fa')
    output:
        ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz','.upper.fa')
    shell:
        "awk 'BEGIN{{FS=\" \"}}{{if(!/>/){{print toupper($0)}}else{{print $1}}}}' {input} > {output}"

rule index_reference_fasta:
    input:
        ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa")
    output:
        ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa.fai")
    conda:
        'envs/samtools.yaml'
    resources:
        runtime=62,
        mem_mb = 8 * 1024
    shell:
        "samtools faidx {input}"

rule download_par_regions:
    output:
        ref_dir + par_regions_bed_url.split('/')[-1]
    params:
        url=par_regions_bed_url,
        output_dir=ref_dir
    localrule: True
    shell:
        "wget -q -P {params.output_dir} {params.url}"


def get_barcoded_pb_bam_file(wildcards):
    barcoded_bam = glob.glob(f'{data_dir}{wildcards.sample}/hifi_reads/*[!unassigned].bam')
    return barcoded_bam

def get_barcoded_ont_bam_file(wildcards):
    barcoded_bam = glob.glob(f'{data_dir}{wildcards.sample}/ont_r10.4.1_dorado/*.bam')
    return barcoded_bam

def get_barcoded_ul_ont_fq_file(wildcards):
    barcoded_bam = glob.glob(f'{data_dir}{wildcards.sample}/ul_ont_r10.4.1_dorado/*.fastq.gz')
    return barcoded_bam


rule bam2fq_pacbio:
    input:
        get_barcoded_pb_bam_file
    output:
        base_results_dir + "{sample}/reads/hifi_reads.fq.gz"
    conda:
        "envs/samtools.yaml"
    retries: 3
    resources:
        runtime=lambda wildcards, attempt: (attempt+1) * 180,
        mem_mb=16 * 1024
    threads: 4
    shell:
        "for bam in {input}; do samtools fastq -@ {threads} -T MM,ML ${{bam}} | bgzip -c >> {output}; done"

rule bam2fq_ont:
    input:
        get_barcoded_ont_bam_file
    output:
        base_results_dir + "{sample}/reads/ont_r10.4.1_dorado_reads.fq.gz"
    conda:
        "envs/samtools.yaml"
    resources:
        runtime=480,
        mem_mb=32 * 1024
    threads: 16
    shell:
        "for bam in {input}; do samtools fastq -@ {threads} -T MM,ML ${{bam}} | bgzip -c >> {output}; done"


rule merge_fq_ul_ont:
    input:
        get_barcoded_ul_ont_fq_file
    output:
        base_results_dir + "{sample}/reads/ul_ont_r10.4.1_dorado_reads.fq.gz"
    conda:
        "envs/samtools.yaml"
    resources:
        runtime=480,
        mem_mb=32 * 1024
    threads: 16
    shell:
        "for fq in {input}; do cat ${{fq}} >> {output}; done"

rule index_fastq:
    input:
        base_results_dir + "{sample}/reads/{seq_tech}_reads.fq.gz"
    output:
        base_results_dir + "{sample}/reads/{seq_tech}_reads.fq.gz.fai"
    conda:
        "envs/samtools.yaml"
    resources:
        runtime=2 * 60,
        mem_mb=16 * 1024
    threads: 1
    shell:
        "samtools faidx {input}"

def get_input_align_reads(wildcards):
    if wildcards.aligner == 'minimap2':
        return {"reads": base_results_dir + f"{wildcards.sample}/reads/{wildcards.seq_tech}_reads.fq.gz",
                "ref": ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa"),
                "fai": ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa.fai")}
    elif wildcards.aligner == 'winnowmap':
        rep_kmer=base_results_dir + f"{wildcards.sample}/alignments/{wildcards.assembler}.{wildcards.asm}.repetitive_k15.txt"
        return {"reads": base_results_dir + f"{wildcards.sample}/reads/{wildcards.seq_tech}_reads.fq.gz",
                "ref": ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa"),
                "fai": ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa.fai"),
                "rep_kmer": rep_kmer}
    elif wildcards.aligner == 'pbmm2':
        return {"reads": glob.glob(f'{data_dir}{wildcards.sample}/{wildcards.seq_tech}_reads/*[!unassigned].bam')[0],
                "ref": ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa"),
                "fai": ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa.fai")}


def get_aligner_command(wildcards):
    if wildcards.aligner == 'minimap2' and wildcards.seq_tech == 'hifi':
        return "minimap2 -ax map-hifi -I8g -p0.5 --eqx --cs -Y -L -k 19 --MD -y -t"
    elif wildcards.aligner == 'minimap2' and 'ont' in wildcards.seq_tech:
        return "minimap2 -ax map-ont -I8g -p0.5 --eqx --cs -Y -L -k 19 --MD -y -t"
    elif wildcards.aligner == 'winnowmap' and wildcards.seq_tech == 'hifi':
        return "winnowmap -W {input.rep_kmer} -ax map-pb -I8g -p0.5 --eqx --cs -Y -L --MD -y -t"
    elif wildcards.aligner == 'winnowmap' and 'ont' in wildcards.seq_tech:
        return "winnowmap -W {input.rep_kmer} -ax map-ont -I8g -p0.5 --eqx --cs -Y -L --MD -y -t"
    elif wildcards.aligner == "pbmm2" and wildcards.seq_tech == 'hifi':
        return "pbmm2 align --preset HiFi -j"

rule align_reads_to_ref:
    input:
        unpack(get_input_align_reads)
    output:
        bam = base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.bam",
        bai= base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.bam.bai"
    params:
        temp_dir = base_results_dir + "{sample}/alignments/",
        aligner_cmd = get_aligner_command
    conda:
        "envs/pbmm2.yaml"
    threads: 32
    resources:
        runtime=24 * 60,
        mem_mb = 64 * 1024
    wildcard_constraints:
        aligner='pbmm2|minimap2|winnowmap',
        seq_tech= '|'.join(seq_tech)
    shell:
        "{params.aligner_cmd} {threads} {input.ref} {input.reads} | "
        "samtools sort -@ {threads} -T {params.temp_dir} | "
        "samtools view -b -@ {threads} -T {params.temp_dir} > {output.bam};"
        "samtools index -b -@ {threads} {output.bam}"

rule calculate_depth_on_Y_chrom:
    input:
        base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.bam"
    output:
        base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.depth_Y_chrom.txt"
    conda:
        "envs/samtools.yaml"
    resources:
        runtime=4 * 60,
        mem_mb=8 * 1024
    threads: 4
    shell:
        "samtools depth -@ {threads} -aa -r chrY {input} | awk '{{sum += $3}}END{{print sum / NR}}' > {output}"

rule install_deepvariant:
    output:
        "docker/deepvariant_latest.sif"
    localrule: True
    shell:
        "singularity pull docker://google/deepvariant; mkdir -p docker/; mv deepvariant_latest.sif docker/"

def get_deepvariant_cmd(wildcards, input, output):
    with open(input.depth_y,'r') as f:
        depth = float(f.readline().strip())
    f.close()
    if wildcards.chrom != 'Y' and wildcards.chrom != 'X':
        return f"singularity exec --bind {base_results_dir} --bind /usr/lib/locale/ {input.sif} "+\
               f"/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref {input.ref} --reads {input.bam} "+\
               f"--regions chr{wildcards.chrom} --output_vcf {output.vcf} --num_shards 32 "+\
               f"--intermediate_results_dir {output.tmp}"
    # male
    if depth >= approx_sequenced_coverage / 4:
        if wildcards.chrom == 'Y' or wildcards.chrom == 'X':
            return f"singularity exec --bind {base_results_dir} --bind /usr/lib/locale/ {input.sif} " + \
                   f"/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref {input.ref} --reads {input.bam} " + \
                   f"--regions chr{wildcards.chrom} --output_vcf {output.vcf} --num_shards 32 " + \
                   f"--intermediate_results_dir {output.tmp} --haploid_contigs \"chrX,chrY\" --par_regions_bed {input.par}"
    # female
    else:
         if wildcards.chrom == "X":
             return f"singularity exec --bind {base_results_dir} --bind /usr/lib/locale/ {input.sif} " + \
                    f"/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref {input.ref} --reads {input.bam} " + \
                    f"--regions chr{wildcards.chrom} --output_vcf {output.vcf} --num_shards 32 --intermediate_results_dir {output.tmp}"
         # Y chromosome does not exist
         else:
             return f"touch {output.vcf}; touch {output.tbi}; mkdir -p {output.tmp}"


rule run_deepvariant:
    input:
        bam=base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.bam",
        ref=ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa"),
        fai = ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa.fai"),
        sif= "docker/deepvariant_latest.sif",
        par= ref_dir + par_regions_bed_url.split('/')[-1],
        depth_y=base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.depth_Y_chrom.txt"
    output:
        vcf=base_results_dir + "{sample}/deepvariant/deepvariant_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz",
        tbi=base_results_dir + "{sample}/deepvariant/deepvariant_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz.tbi",
        tmp=temp(directory(base_results_dir + "{sample}/deepvariant/tmp_{seq_tech}_deepvariant_{aligner}_{chrom}"))
    params:
        cmd=lambda wildcards, input, output: get_deepvariant_cmd(wildcards, input, output)
    conda:
        "envs/samtools.yaml"
    wildcard_constraints:
        aligner=aligner,
        chrom="|".join(autosomes + sex_chromosomes)
    threads: 32
    resources:
        mem_mb=48 * 1024,
        runtime=24 * 60
    shell:
        "{params.cmd}"

rule filter_variants:
    input:
        base_results_dir + "{sample}/deepvariant/deepvariant_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz"
    output:
        vcf=base_results_dir + "{sample}/deepvariant/deepvariant_{seq_tech}_hs1.{aligner}.chr{chrom}.filtered.vcf.gz",
        tbi=base_results_dir + "{sample}/deepvariant/deepvariant_{seq_tech}_hs1.{aligner}.chr{chrom}.filtered.vcf.gz.tbi"
    resources:
        runtime=30,
        mem_mb=4 * 1024
    conda:
        "envs/samtools.yaml"
    wildcard_constraints:
        aligner=aligner
    shell:
        "bcftools view -f \"PASS\" -i \"GQ>20\" -Oz -o {output.vcf} {input}; "
        "tabix {output.vcf}"


rule run_sniffles:
    input:
        bam=base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.bam",
        ref=ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa")
    output:
        vcf=temp(base_results_dir + "{sample}/sniffles/sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}_tmp.vcf.gz")
    threads: 16
    resources:
        runtime=240,
        mem_mb=32 * 1024
    params:
        chrom="chr{chrom}"
    wildcard_constraints:
        chrom="|".join(autosomes + sex_chromosomes)
    conda:
        "envs/sniffles.yaml"
    shell:
        "sniffles -i {input.bam} --vcf {output.vcf} --reference {input.ref} -t {threads} -c {params.chrom} --output-rnames"

rule reheader_sniffls:
    input:
        base_results_dir + "{sample}/sniffles/sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}_tmp.vcf.gz"
    output:
        vcf=base_results_dir + "{sample}/sniffles/sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz",
        header=temp(base_results_dir + "{sample}/sniffles/sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.header")
    localrule: True
    wildcard_constraints:
        chrom="|".join(autosomes + sex_chromosomes)
    shell:
        "bcftools view -h {input} | sed 's/SAMPLE/default/' > {output.header}; "
        "bcftools reheader -h {output.header} {input} > {output.vcf}; tabix {output.vcf}"


rule filter_sniffles:
    input:
        base_results_dir + "{sample}/sniffles/sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz"
    output:
        base_results_dir + "{sample}/sniffles/sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.filtered.vcf.gz"
    conda:
        "envs/samtools.yaml"
    resources:
        runtime=10,
        mem_mb=4 * 1024
    wildcard_constraints:
        chrom="|".join(autosomes + sex_chromosomes)
    shell:
        "bcftools view -i 'SVTYPE!=\"BND\" && GQ >=20' -Oz -o {output} {input}; tabix {output}"

rule merge_snvs_svs:
    input:
        dv=base_results_dir + "{sample}/deepvariant/deepvariant_{seq_tech}_hs1.{aligner}.chr{chrom}.filtered.vcf.gz",
        dv_tbi=base_results_dir + "{sample}/deepvariant/deepvariant_{seq_tech}_hs1.{aligner}.chr{chrom}.filtered.vcf.gz.tbi",
        sniffles=base_results_dir + "{sample}/sniffles/sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.filtered.vcf.gz",
        ref=ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa")
    output:
        vcf = base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz",
        tbi = base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz.tbi"
    conda:
        "envs/samtools.yaml"
    retries: 3
    resources:
        runtime=lambda wildcards, attempt: (attempt + 1) * 62,
        mem_mb=4 * 1024
    wildcard_constraints:
        chrom="|".join(autosomes + sex_chromosomes)
    shell:
        "bcftools concat {input.dv} {input.sniffles} | "
        "bcftools sort | "
        "bcftools norm -m+ -f {input.ref} |"
        "bcftools sort -Oz -o {output.vcf} ; tabix {output.vcf}"


rule call_methylations:
    input:
        bam=base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.bam",
    output:
        base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.methylation.combined.bed.gz"
    params:
        prefix=base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.methylation",
        exec=aligned_bam_to_cpg_scores_path
    threads: 32
    resources:
        runtime=120,
        mem_mb=32 * 1024
    wildcard_constraints:
        aligner=aligner
    shell:
        "{params.exec} --bam {input.bam} --output-prefix {params.prefix} --threads {threads}"


def get_longhap_options(wildcards):
    if 'ont' in wildcards.seq_tech:
        return "--ont"
    else:
        return "--pacbio"

rule run_longhap_meth_phase_indels:
    input:
        vcf=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz",
        tbi=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz.tbi",
        bam=base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.bam",
        ref=ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa"),
        methylation=base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.methylation.combined.bed.gz"
    output:
        vcf=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_phased.vcf.gz",
        blocks=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_blocks.bed",
        npz=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_transition_mat_variants.npz",
        npz_meth=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_transition_mat_meth.npz",
        read_states=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_read_states.json",
        meth=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_diff_meth_sites.tab",
        variant_read_mapping=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_variant_read_mapping.json",
        unphaseable_variants=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_unphaseable_variants.npz",
        allele_coverage=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth_allele_coverage.npz"
    conda:
        "envs/longhap.yaml"
    params:
        chrom="chr{chrom}",
        options=get_longhap_options,
        longhap_path=longhap_path
    wildcard_constraints:
        aligner=aligner,
        seq_tech="|".join(seq_tech),
        chrom="|".join(autosomes + sex_chromosomes)
    threads: 1
    resources:
        runtime=30,
        mem_mb=16 * 1024
    benchmark: base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_meth.benchmark"
    shell:
        "{params.longhap_path} --vcf {input.vcf} -b {input.bam} -c {params.chrom} -r {input.ref} "
        "-o {output.vcf} --output_blocks {output.blocks} --output_transition_matrix {output.npz} "
        "--output_variant_read_mapping {output.variant_read_mapping} "
        "--output_unphaseable_variants {output.unphaseable_variants} "
        "--output_read_states {output.read_states} --output_allele_coverage {output.allele_coverage} " 
        "-m {input.methylation} --output_transition_matrix_meth {output.npz_meth} "
        "--output_differentially_methylated_sites {output.meth} {params.options} --force"

rule run_longhap_phase_indels:
    input:
        vcf=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz",
        tbi=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.vcf.gz.tbi",
        bam=base_results_dir + "{sample}/alignments/{seq_tech}_reads_to_hs1.{aligner}.bam",
        ref=ref_dir + reference_genome_url.split('/')[-1].replace('.fa.gz', ".upper.fa")
    output:
        vcf=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_phased.vcf.gz",
        blocks=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_blocks.bed",
        npz=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_transition_mat.npz",
        read_states=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_read_states.json",
        variant_read_mapping=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_variant_read_mapping.json",
        unphaseable_variants=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_unphaseable_variants.npz",
        allele_coverage=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap_allele_coverage.npz"
    conda:
        "envs/longhap.yaml"
    params:
        chrom="chr{chrom}",
        options=get_longhap_options,
        longhap_path=longhap_path
    wildcard_constraints:
        aligner=aligner,
        seq_tech="|".join(seq_tech),
        chrom="|".join(autosomes + sex_chromosomes)
    threads: 1
    resources:
        runtime=30,
        mem_mb=16 * 1024
    benchmark: base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.longhap.benchmark"
    shell:
        "{params.longhap_path} --vcf {input.vcf} -b {input.bam} -c {params.chrom} -r {input.ref} "
        "-o {output.vcf} --output_blocks {output.blocks} --output_transition_matrix {output.npz} "
        "--output_read_states {output.read_states} --output_allele_coverage {output.allele_coverage} " 
        "--output_variant_read_mapping {output.variant_read_mapping} "
        "--output_unphaseable_variants {output.unphaseable_variants} {params.options} --force"

rule index_phased_vcf:
    input:
        vcf=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.vcf.gz"
    output:
        base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.vcf.gz.tbi"
    conda:
        "envs/samtools.yaml"
    resources:
        runtime=10,
        mem_mb=2 * 1024
    wildcard_constraints:
        aligner=aligner,
        seq_tech="|".join(seq_tech),
        chrom="|".join(autosomes + sex_chromosomes)
    shell:
        "tabix {input}"


rule run_whatshap_stats:
    input:
        base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.vcf.gz"
    output:
        gtf=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.stats.gtf",
        tsv=base_results_dir+ "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.stats.tsv"
    conda:
        "envs/whatshap.yaml"
    resources:
        runtime=20,
        mem_mb=4 * 1024
    wildcard_constraints:
        aligner=aligner,
        seq_tech="|".join(seq_tech)
    shell:
        'whatshap stats --gtf={output.gtf} --tsv={output.tsv} {input}'

rule run_whatshap_stats_complex:
    input:
        base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.vcf.gz"
    output:
        gtf=base_results_dir + "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.stats_complex_variants.gtf",
        tsv=base_results_dir+ "{sample}/longhap/deepvariant_sniffles_{seq_tech}_hs1.{aligner}.chr{chrom}.{phasing_method}_phased.stats_complex_variants.tsv"
    conda:
        "envs/whatshap.yaml"
    resources:
        runtime=20,
        mem_mb=4 * 1024
    wildcard_constraints:
        aligner=aligner,
        seq_tech="|".join(seq_tech)
    shell:
        'whatshap stats --gtf={output.gtf} --tsv={output.tsv} <(bcftools view -V snps {input})'
