configfile: "config.yaml"

from glob import glob
import os

FASTQ_DIR = config["samples_dir"]

SAMPLES = sorted({
    os.path.basename(f).split("_R")[0]
    for f in glob(f"{FASTQ_DIR}/*_R1.fastq.gz")
})

rule all:
    input:
        "results/multiqc/multiqc_report.html",
        "results/counts/gene_counts.txt"

rule fastqc_raw:
    input:
        fq=lambda wildcards: f"{FASTQ_DIR}/{wildcards.sample}_{wildcards.read}.fastq.gz"
    output:
        zip="results/qc/fastqc_raw/{sample}_{read}_fastqc.zip",
        html="results/qc/fastqc_raw/{sample}_{read}_fastqc.html"
    log:
        "logs/fastqc_raw/{sample}_{read}.log"
    threads: config["threads"]["fastqc"]
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        mkdir -p results/qc/fastqc_raw logs/fastqc_raw
        fastqc {input.fq} \
            --outdir results/qc/fastqc_raw \
            --threads {threads} \
            > {log} 2>&1
        """

rule trim_galore:
    input:
        r1=lambda wc: f"{FASTQ_DIR}/{wc.sample}_R1.fastq.gz",
        r2=lambda wc: f"{FASTQ_DIR}/{wc.sample}_R2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_R1_val_1.fq.gz",
        r2="results/trimmed/{sample}_R2_val_2.fq.gz"
    log:
        "logs/trim_galore/{sample}.log"
    threads: config["threads"]["trim"]
    conda:
        "envs/trim.yaml"
    shell:
        """
        mkdir -p results/trimmed logs/trim_galore
        trim_galore \
            --paired \
            --cores {threads} \
            --output_dir results/trimmed \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """

rule fastqc_trimmed:
    input:
        fq=lambda wc: f"results/trimmed/{wc.sample}_{wc.read}_val_{1 if wc.read=='R1' else 2}.fq.gz"
    output:
        directory("results/qc/fastqc_trimmed/{sample}_{read}")
    log:
        "logs/fastqc/trimmed/{sample}_{read}.log"
    threads: config["threads"]["fastqc"]
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        mkdir -p results/qc/fastqc_trimmed logs/fastqc/trimmed
        fastqc {input.fq} \
            --outdir results/qc/fastqc_trimmed \
            --threads {threads} \
            > {log} 2>&1
        """

rule multiqc:
    input:
        expand("results/qc/fastqc_raw/{sample}_{read}_fastqc.zip",
               sample=SAMPLES, read=["R1", "R2"]),
        expand("results/qc/fastqc_trimmed/{sample}_{read}",
               sample=SAMPLES, read=["R1", "R2"])
    output:
        "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        mkdir -p results/multiqc
        multiqc results/qc \
            --outdir results/multiqc \
            > {log} 2>&1
        """

rule star_index:
    input:
        fasta=config["genome_fasta"],
        gtf=config["gtf"]
    output:
        directory(config["star_index"])
    threads: config["threads"]["star"]
    conda:
        "envs/star.yaml"
    shell:
        """
        mkdir -p {output}
        STAR \
          --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {output} \
          --genomeFastaFiles {input.fasta} \
          --sjdbGTFfile {input.gtf}
        """

rule star_align:
    input:
        index=config["star_index"],
        r1="results/trimmed/{sample}_R1_val_1.fq.gz",
        r2="results/trimmed/{sample}_R2_val_2.fq.gz"
    output:
        bam="results/aligned/{sample}.sorted.bam",
        bai="results/aligned/{sample}.sorted.bam.bai"
    log:
        "logs/star/{sample}.log"
    threads: config["threads"]["star"]
    conda:
        "envs/star.yaml"
    shell:
        """
        mkdir -p results/aligned logs/star
        STAR \
          --genomeDir {input.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --runThreadN {threads} \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix results/aligned/{wildcards.sample}. \
          > {log} 2>&1

        mv results/aligned/{wildcards.sample}.Aligned.sortedByCoord.out.bam {output.bam}
        samtools index {output.bam}
        """

rule featurecounts:
    input:
        bam=expand(
            "results/aligned/{sample}.sorted.bam",
            sample=SAMPLES
        ),
        gtf=config["gtf"]
    output:
        counts="results/counts/gene_counts.txt",
        summary="results/counts/gene_counts.txt.summary"
    log:
        "logs/featurecounts.log"
    threads: config["threads"]["counts"]
    conda:
        "envs/featurecounts.yaml"
    params:
        paired="-p" if config["library"]["paired"] else "",
        strand=config["library"]["stranded"],
        feature=config["featurecounts"]["feature"],
        attribute=config["featurecounts"]["attribute"]

    shell:
        """
        mkdir -p results/counts logs
        featureCounts \
            -T {threads} \
            {params.paired} \
            -s {params.strand} \
            -t {params.feature} \
            -g {params.attribute} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bam} \
            > {log} 2>&1
        """
