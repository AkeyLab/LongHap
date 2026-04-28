import numpy as np
base_dir = '/scratch/gpfs/AKEY/apfennig/LongHap/'
workdir: base_dir

longhap_vcf =  base_dir + "{sample}/longhap/deepvariant_sniffles_hifi_hs1.minimap2.chr{chrom}.{tool}_phased.vcf.gz"
reference_vcf =  "/scratch/gpfs/AKEY/apfennig/1KGP_T2T_Rhie_et_al_Nature_2023/chr{chrom}.recalibrated.snp_indel.pass.vcf.gz"

maternal_assembly = base_dir + "{sample}/assemblies/maternal_assembly.fa.gz"
paternal_assembly = base_dir + "{sample}/assemblies/paternal_assembly.fa.gz"

reference_fasta = base_dir + 'reference/hs1.upper.fa'
reference_fasta_index = base_dir + 'reference/hs1.upper.fa.fai'
par=base_dir + 'reference/chm13v2.0_PAR.bed'

omim_regions=base_dir + "reference/genemap2_t2t.txt"

samples = ['HG00609', 'HG00658', 'HG00738', 'HG01099', 'HG02723', 'HG02615', 'HG002']
samples_sex = {'HG00609': "male",
               'HG00658': "male",
               'HG00738': "male",
               'HG01099': "male",
               'HG02723': "female",
               'HG02615': "female",
               'HG002': "male"}

chromosomes = np.arange(1, 23)

tools = ['longhap', 'longhap_meth']

rule all:
    input:
        expand(longhap_vcf.replace('.vcf.gz','.annotated.stats.tsv'),sample=samples,
            chrom=chromosomes,tool=tools),
        expand(longhap_vcf.replace('.vcf.gz', '.annotated.rare.stats.tsv'), sample=samples,
            chrom=chromosomes, tool=tools),
        expand(longhap_vcf.replace('.vcf.gz','.annotated.complex.stats.tsv'), sample=samples,
            chrom=chromosomes, tool=tools),
        expand(longhap_vcf.replace('.vcf.gz', '.annotated.rare.complex.stats.tsv'), sample=samples,
            chrom=chromosomes, tool=tools),
        expand(longhap_vcf.replace('.vcf.gz', '.annotated.evaluation.tab'), sample=samples,
            chrom=chromosomes, tool=tools),
        expand(longhap_vcf.replace('.vcf.gz', '.annotated.complex.evaluation.tab'), sample=samples,
            chrom=chromosomes, tool=tools),
        expand(longhap_vcf.replace('.vcf.gz','.annotated.rare.evaluation.tab'),sample=samples,
            chrom=chromosomes, tool=tools),
        expand(longhap_vcf.replace('.vcf.gz','.annotated.rare.complex.evaluation.tab'),sample=samples,
            chrom=chromosomes, tool=tools),
        expand(base_dir + "{tool}_switch_errors_per_10kb_rare_complex.bed", tool=tools),
        expand(base_dir + "{tool}_switch_errors_per_10kb_complex.bed", tool=tools),
        expand(base_dir + "{tool}_switch_errors_per_10kb_rare.bed", tool=tools),
        expand(base_dir + "{tool}_switch_errors_per_10kb.bed", tool=tools),
        expand(base_dir + "{sample}/longhap/{tool}_omim_regions_switch_error_counts.bed", sample=samples, tool=tools),
        expand(base_dir + "{sample}/longhap/{tool}_omim_regions_switch_error_counts_rare.bed", sample=samples, tool=tools),
        expand(base_dir + "{sample}/longhap/{tool}_omim_regions_switch_error_counts_complex.bed", sample=samples, tool=tools),
        expand(base_dir + "{sample}/longhap/{tool}_omim_regions_switch_error_counts_rare_complex.bed", sample=samples, tool=tools)

rule window_reference:
    input:
        reference_fasta_index
    output:
        reference_fasta_index.replace('.fa.fai', '.10kb_windows.bed')
    localrule: True
    shell:
        "bedtools makewindows -g <(cut -f1,2 {input}) -w 10000 -s 10000 > {output}"

rule norm_reference_vcf:
    input:
        reference_vcf
    output:
        reference_vcf.replace('.vcf.gz', '.normalized.vcf.gz')
    resources:
        mem_mb=8000,
        runtime=18*60
    conda:
        "envs/samtools.yaml"
    threads: 8
    shell:
        "bcftools norm --threads {threads} -m -any {input} -Oz -o {output}; "
        "tabix {output}"

rule norm_longhap_vcf:
    input:
        longhap_vcf
    output:
        longhap_vcf.replace('.vcf.gz', '.normalized.vcf.gz')
    resources:
        mem_mb=4000,
        runtime=120
    conda:
        "envs/samtools.yaml"
    threads: 4
    shell:
        "bcftools norm --threads {threads} -m -any {input} -Oz -o {output}; "
        "tabix {output}"

rule annotate_vcf_with_allele_frequencies:
    input:
        phased_vcf = longhap_vcf.replace('.vcf.gz', '.normalized.vcf.gz'),
        reference_vcf = reference_vcf.replace('.vcf.gz', '.normalized.vcf.gz')
    output:
        annotated_vcf = longhap_vcf.replace('.vcf.gz', '.annotated.vcf.gz')
    resources:
        mem_mb=4000,
        runtime=480
    conda:
        "envs/samtools.yaml"
    threads: 4
    shell:
        "bcftools annotate --threads {threads} -a {input.reference_vcf} -c INFO/AF "
        "-Oz -o {output.annotated_vcf} {input.phased_vcf}; "
        "tabix {output.annotated_vcf}"

rule get_rare_variants:
    input:
        longhap_vcf.replace('.vcf.gz','.annotated.vcf.gz')
    output:
        rare_vcf = longhap_vcf.replace('.vcf.gz', '.annotated.rare.vcf.gz')
    resources:
        mem_mb=4000,
        runtime=120
    conda:
        "envs/samtools.yaml"
    threads: 4
    shell:
        "bcftools view --threads {threads} -i 'INFO/AF<0.01' -Oz -o {output.rare_vcf} {input}; "
        "tabix {output.rare_vcf}"


rule run_whatshap_stats:
    input:
        rare_vcf = longhap_vcf.replace('.vcf.gz', '.annotated.vcf.gz')
    output:
        gtf=longhap_vcf.replace('.vcf.gz', '.annotated.stats.gtf'),
        tsv=longhap_vcf.replace('.vcf.gz', '.annotated.stats.tsv')
    conda:
        "envs/whatshap.yaml"
    resources:
        runtime=20,
        mem_mb=4 * 1024
    shell:
        'whatshap stats --gtf={output.gtf} --tsv={output.tsv} {input}'

rule run_whatshap_stats_rare:
    input:
        rare_vcf = longhap_vcf.replace('.vcf.gz', '.annotated.rare.vcf.gz')
    output:
        gtf=longhap_vcf.replace('.vcf.gz', '.annotated.rare.stats.gtf'),
        tsv=longhap_vcf.replace('.vcf.gz', '.annotated.rare.stats.tsv')
    conda:
        "envs/whatshap.yaml"
    resources:
        runtime=20,
        mem_mb=4 * 1024
    shell:
        'whatshap stats --gtf={output.gtf} --tsv={output.tsv} {input}'

rule run_whatshap_stats_complex:
    input:
        rare_vcf = longhap_vcf.replace('.vcf.gz', '.annotated.vcf.gz')
    output:
        gtf = longhap_vcf.replace('.vcf.gz','.annotated.complex.stats.gtf'),
        tsv = longhap_vcf.replace('.vcf.gz','.annotated.complex.stats.tsv')
    conda:
        "envs/whatshap.yaml"
    resources:
        runtime=20,
        mem_mb=4 * 1024
    shell:
        'whatshap stats --gtf={output.gtf} --tsv={output.tsv} <(bcftools view -V snps {input})'


rule run_whatshap_stats_rare_complex:
    input:
        rare_vcf = longhap_vcf.replace('.vcf.gz', '.annotated.rare.vcf.gz')
    output:
        gtf = longhap_vcf.replace('.vcf.gz','.annotated.rare.complex.stats.gtf'),
        tsv = longhap_vcf.replace('.vcf.gz','.annotated.rare.complex.stats.tsv')
    conda:
        "envs/whatshap.yaml"
    resources:
        runtime=20,
        mem_mb=4 * 1024
    shell:
        'whatshap stats --gtf={output.gtf} --tsv={output.tsv} <(bcftools view -V snps {input})'

def get_sex_specific_dipcall_options(wildcards):
    options = ''
    if samples_sex[wildcards.sample] == 'male':
        options += f'-x {par}'
    return options

rule call_dipcall:
    input:
        pat=paternal_assembly,
        mat =maternal_assembly,
        ref = reference_fasta,
        ref_idx=reference_fasta_index
    output:
        base_dir + "{sample}/dipcall/asm.dip.vcf.gz"
    params:
        prefix=base_dir + "{sample}/dipcall/asm",
        output_dir =base_dir + "{sample}/dipcall/",
        sex_specific_options = get_sex_specific_dipcall_options
    threads: 16
    resources:
        mem_mb = 128 * 1024,
        runtime=10 * 60
    shell:
        "mkdir -p {params.output_dir}; "
        "rm -r {params.output_dir}; "
        "mkdir -p {params.output_dir}; "
        "run-dipcall -t {threads} {params.sex_specific_options} {params.prefix} {input.ref} {input.pat} {input.mat} > {params.prefix}.mak;"
        "make -j{threads} -f {params.prefix}.mak"

rule index_assembly_calls:
    input:
        base_dir + "{sample}/dipcall/asm.dip.vcf.gz"
    output:
        base_dir + "{sample}/dipcall/asm.dip.vcf.gz.tbi"
    conda:
        "envs/samtools.yaml"
    resources:
        mem_mb=4000,
        runtime=62
    threads: 1
    shell:
        "tabix {input}"

rule get_phased_sites:
    input:
        vcf=longhap_vcf.replace('.vcf.gz', '.annotated.vcf.gz'),
        tbi=longhap_vcf.replace('.vcf.gz', '.annotated.vcf.gz')
    output:
        tab=longhap_vcf.replace('.vcf.gz', '.annotated.phased_sites.tab')
    conda:
        "envs/samtools.yaml"
    resources:
        mem_mb=4 * 1024,
        runtime=20
    shell:
        "bcftools view -g het -p -H {input.vcf} | cut -f1,2 > {output.tab}"

rule get_phased_rare_sites:
    input:
        vcf=longhap_vcf.replace('.vcf.gz', '.annotated.rare.vcf.gz'),
        tbi=longhap_vcf.replace('.vcf.gz', '.annotated.rare,vcf.gz')
    output:
        tab=longhap_vcf.replace('.vcf.gz', '.annotated.rare.phased_sites.tab')
    conda:
        "envs/samtools.yaml"
    resources:
        mem_mb=4 * 1024,
        runtime=20
    shell:
        "bcftools view -g het -p -H {input.vcf} | cut -f1,2 > {output.tab}"


def get_input_shared_phased_sites(wildcards):
    return [longhap_vcf.replace('.vcf.gz', '.annotated.phased_sites.tab').format(sample=wildcards.sample,
        chrom=wildcards.chrom, tool=tool) for tool in tools]

rule get_shared_phased_sites:
    input:
        get_input_shared_phased_sites
    output:
        base_dir + "{sample}/longhap/shared_phased_sites_chr{chrom}.tab"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    params:
        n_phasing_methods = lambda wildcards: len(tools)
    shell:
        "cat {input} | sort | uniq -c | "
        "awk 'BEGIN{{OFS=\"\\t\"}}{{if ($1 == {params.n_phasing_methods}) print $2, $3}}' | sort -k2n > {output}"



def get_input_shared_phased_rare_sites(wildcards):
    return [longhap_vcf.replace('.vcf.gz', '.annotated.rare.phased_sites.tab').format(sample=wildcards.sample,
        chrom=wildcards.chrom, tool=tool) for tool in tools]

rule get_shared_phased_rare_sites:
    input:
        get_input_shared_phased_sites
    output:
        base_dir + "{sample}/longhap/shared_phased_rare_sites_chr{chrom}.tab"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    params:
        n_phasing_methods = lambda wildcards: len(tools)
    shell:
        "cat {input} | sort | uniq -c | "
        "awk 'BEGIN{{OFS=\"\\t\"}}{{if ($1 == {params.n_phasing_methods}) print $2, $3}}' | sort -k2n > {output}"


rule run_whatshap_compare:
    input:
        vcf=longhap_vcf.replace('.vcf.gz', '.annotated.vcf.gz'),
        gt=base_dir + "{sample}/dipcall/asm.dip.vcf.gz",
        tbi=base_dir + "{sample}/dipcall/asm.dip.vcf.gz.tbi",
        shared_sites=base_dir + "{sample}/longhap/shared_phased_sites_chr{chrom}.tab"
    output:
        bed=longhap_vcf.replace('.vcf.gz', '.annotated.switch_errors.bed'),
        tab=longhap_vcf.replace('.vcf.gz', '.annotated.evaluation.tab')
    conda:
        "envs/whatshap.yaml"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    params:
        chrom='chr{chrom}'
    shell:
        "whatshap compare --names truth,query --ignore-sample-name --tsv-pairwise {output.tab} "
        "--switch-error-bed {output.bed} <(bcftools norm -r {params.chrom} -m -any {input.gt}) "
        "<(bcftools annotate -x \"FORMAT/PS\" {input.vcf} | bcftools norm -m -any | "
        "bcftools view -T <(sort -k2n {input.shared_sites}))"


rule run_whatshap_compare_rare:
    input:
        vcf=longhap_vcf.replace('.vcf.gz', '.annotated.rare.vcf.gz'),
        gt=base_dir + "{sample}/dipcall/asm.dip.vcf.gz",
        tbi=base_dir + "{sample}/dipcall/asm.dip.vcf.gz.tbi",
        shared_sites=base_dir + "{sample}/longhap/shared_phased_rare_sites_chr{chrom}.tab"
    output:
        bed=longhap_vcf.replace('.vcf.gz', '.annotated.rare.switch_errors.bed'),
        tab=longhap_vcf.replace('.vcf.gz', '.annotated.rare.evaluation.tab')
    conda:
        "envs/whatshap.yaml"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    params:
        chrom='chr{chrom}'
    shell:
        "whatshap compare --names truth,query --ignore-sample-name --tsv-pairwise {output.tab} "
        "--switch-error-bed {output.bed} <(bcftools norm -r {params.chrom} -m -any {input.gt}) "
        "<(bcftools annotate -x \"FORMAT/PS\" {input.vcf} | bcftools norm -m -any | "
        "bcftools view -T <(sort -k2n {input.shared_sites}))"


rule run_whatshap_compare_complex:
    input:
        vcf=longhap_vcf.replace('.vcf.gz', '.annotated.vcf.gz'),
        gt=base_dir + "{sample}/dipcall/asm.dip.vcf.gz",
        tbi=base_dir+ "{sample}/dipcall/asm.dip.vcf.gz.tbi",
        shared_sites=base_dir + "{sample}/longhap/shared_phased_sites_chr{chrom}.tab"
    output:
        bed=longhap_vcf.replace('.vcf.gz', '.annotated.complex.switch_errors.bed'),
        tab=longhap_vcf.replace('.vcf.gz', '.annotated.complex.evaluation.tab')
    conda:
        "envs/whatshap.yaml"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    params:
        chrom='chr{chrom}'
    shell:
        "whatshap compare --names truth,query --ignore-sample-name --tsv-pairwise {output.tab} "
        "--switch-error-bed {output.bed} <(bcftools norm -r {params.chrom} -m -any {input.gt}) "
        "<(bcftools annotate -x \"FORMAT/PS\" {input.vcf} | bcftools norm -m -any  | "
        "bcftools view -T <(sort -k2n {input.shared_sites}) | bcftools view -V snps)"


rule run_whatshap_compare_rare_complex:
    input:
        vcf=longhap_vcf.replace('.vcf.gz', '.annotated.rare.vcf.gz'),
        gt=base_dir + "{sample}/dipcall/asm.dip.vcf.gz",
        tbi=base_dir + "{sample}/dipcall/asm.dip.vcf.gz.tbi",
        shared_sites=base_dir + "{sample}/longhap/shared_phased_rare_sites_chr{chrom}.tab"
    output:
        bed=longhap_vcf.replace('.vcf.gz', '.annotated.rare.complex.switch_errors.bed'),
        tab=longhap_vcf.replace('.vcf.gz', '.annotated.rare.complex.evaluation.tab')
    conda:
        "envs/whatshap.yaml"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    params:
        chrom='chr{chrom}'
    shell:
        "whatshap compare --names truth,query --ignore-sample-name --tsv-pairwise {output.tab} "
        "--switch-error-bed {output.bed} <(bcftools norm -r {params.chrom} -m -any {input.gt}) "
        "<(bcftools annotate -x \"FORMAT/PS\" {input.vcf} | bcftools norm -m -any  | "
        "bcftools view -T <(sort -k2n {input.shared_sites}) | bcftools view -V snps)"


def get_switch_errors_omim(wildcards):
    return [longhap_vcf.replace('.vcf.gz','.annotated.switch_errors.bed').format(sample=wildcards.sample,
        tool=wildcards.tool, chrom=chrom) for chrom in chromosomes]

rule omim_region_switch_error_counts:
    input:
        omim_regions=omim_regions,
        switch_errors=get_switch_errors_omim
    output:
        base_dir + "{sample}/longhap/{tool}_omim_regions_switch_error_counts.bed"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    shell:
        "bedtools intersect -a {input.omim_regions} -b <(cat {input.switch_errors} | sort -k1 -k2,2n) -c > {output}"


def get_switch_errors_omim_rare(wildcards):
    return [longhap_vcf.replace('.vcf.gz','.annotated.rare.switch_errors.bed').format(sample=wildcards.sample,
        tool=wildcards.tool, chrom=chrom) for chrom in chromosomes]

rule omim_region_switch_error_counts_rare:
    input:
        omim_regions=omim_regions,
        switch_errors=get_switch_errors_omim_rare
    output:
        base_dir + "{sample}/longhap/{tool}_omim_regions_switch_error_counts_rare.bed"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    shell:
        "bedtools intersect -a {input.omim_regions} -b <(cat {input.switch_errors} | sort -k1 -k2,2n) -c > {output}"

def get_switch_errors_omim_complex(wildcards):
    return [longhap_vcf.replace('.vcf.gz','.annotated.complex.switch_errors.bed').format(sample=wildcards.sample,
        tool=wildcards.tool, chrom=chrom) for chrom in chromosomes]

rule omim_region_switch_error_counts_complex:
    input:
        omim_regions=omim_regions,
        switch_errors=get_switch_errors_omim_complex
    output:
        base_dir + "{sample}/longhap/{tool}_omim_regions_switch_error_counts_complex.bed"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    shell:
        "bedtools intersect -a {input.omim_regions} -b <(cat {input.switch_errors} | sort -k1 -k2,2n) -c > {output}"

def get_switch_errors_omim_rare_complex(wildcards):
    return [longhap_vcf.replace('.vcf.gz','.annotated.rare.complex.switch_errors.bed').format(sample=wildcards.sample,
        tool=wildcards.tool, chrom=chrom) for chrom in chromosomes]

rule omim_region_switch_error_counts_rare_complex:
    input:
        omim_regions=omim_regions,
        switch_errors=get_switch_errors_omim_rare_complex
    output:
        base_dir + "{sample}/longhap/{tool}_omim_regions_switch_error_counts_rare_complex.bed"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    shell:
        "bedtools intersect -a {input.omim_regions} -b <(cat {input.switch_errors} | sort -k1 -k2,2n) -c > {output}"

def get_switch_errors(wildcards):
    return [longhap_vcf.replace('.vcf.gz', '.annotated.switch_errors.bed').format(sample=sample, chrom=chrom,
        tool=wildcards.tool) for sample in samples for chrom in chromosomes]

rule get_switch_errors_per_10kb:
    input:
        windows=reference_fasta_index.replace('.fa.fai', '.10kb_windows.bed'),
        errors=get_switch_errors
    output:
        base_dir + "{tool}_switch_errors_per_10kb.bed"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    shell:
        "bedtools intersect -a {input.windows} -b <({input.errors} | sort -k1 -k2,2n) -c > {output}"


def get_switch_errors_rare(wildcards):
    return [longhap_vcf.replace('.vcf.gz','.annotated.rare.switch_errors.bed').format(sample=sample,chrom=chrom,
        tool=wildcards.tool) for sample in samples for chrom in chromosomes]


rule get_switch_errors_per_10kb_rare:
    input:
        windows=reference_fasta_index.replace('.fa.fai','.10kb_windows.bed'),
        errors=get_switch_errors_rare
    output:
        base_dir + "{tool}_switch_errors_per_10kb_rare.bed"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    shell:
        "bedtools intersect -a {input.windows} -b <({input.errors} | sort -k1 -k2,2n) -c > {output}"


def get_switch_errors_complex(wildcards):
    return [longhap_vcf.replace('.vcf.gz','.annotated.complex.switch_errors.bed').format(sample=sample,chrom=chrom,
        tool=wildcards.tool) for sample in samples for chrom in chromosomes]


rule get_switch_errors_per_10kb_complex:
    input:
        windows=reference_fasta_index.replace('.fa.fai','.10kb_windows.bed'),
        errors=get_switch_errors_complex
    output:
        base_dir + "{tool}_switch_errors_per_10kb_complex.bed"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    shell:
        "bedtools intersect -a {input.windows} -b <({input.errors} | sort -k1 -k2,2n) -c > {output}"

def get_switch_errors_rare_complex(wildcards):
    return [longhap_vcf.replace('.vcf.gz','.annotated.rare.complex.switch_errors.bed').format(sample=sample,chrom=chrom,
        tool=wildcards.tool) for sample in samples for chrom in chromosomes]


rule get_switch_errors_per_10kb_rare_complex:
    input:
        windows=reference_fasta_index.replace('.fa.fai','.10kb_windows.bed'),
        errors=get_switch_errors_rare_complex
    output:
        base_dir + "{tool}_switch_errors_per_10kb_rare_complex.bed"
    resources:
        mem_mb=2 * 1024,
        runtime=62
    shell:
        "bedtools intersect -a {input.windows} -b <({input.errors} | sort -k1 -k2,2n) -c > {output}"


