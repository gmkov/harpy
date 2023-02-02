from snakemake.utils import Paramspace
import pandas as pd

# user specified configs
bam_dir = config["seq_directory"]
contigfile = config["contignames"]
samplenames = config["samplenames"]
variantfile = config["variantfile"]

# declare a dataframe to be a paramspace
paramspace = Paramspace(pd.read_csv(config["paramfile"], sep="\t"), param_sep = "", filename_params="*")

def contignames(contig_file):
    with open(contig_file) as f:
        lines = [line.rstrip() for line in f]
    return lines

contigs = contignames(contigfile)

rule bam_list:
    input: expand(bam_dir + "/{sample}.bam", sample = samplenames)
    output: "Imputation/input/samples.list"
    message: "Creating list of alignment files"
    benchmark: "Benchmark/Impute/filelist.txt"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input:
                fout.write(f"{bamfile}\n")

rule split_contigs:
    input: contigfile
    output: expand("Imputation/input/contigs/{part}", part = contigs)
    message: "Splitting contig names for parallelization"
    benchmark: "Benchmark/Impute/splitcontigs.txt"
    run:
        with open(input[0]) as f:
            cpath = "Imputation/input/contigs"
            for line in f:
                contig = line.rstrip()
                with open(f"{cpath}/{contig}", "w") as fout:
                    gremlin = fout.write(f"{contig}\n")

rule prepare_biallelic_snps:
    input: 
        vcf = variantfile,
        contig = "Imputation/input/contigs/{part}"
    output: pipe("Imputation/input/{part}.bisnp.bcf")
    message: "Keeping only biallelic SNPs from {wildcards.part}"
    benchmark: "Benchmark/Impute/fileprep.{part}.txt"
    threads: 1
    shell:
        """
        bcftools view -m2 -M2 -v snps --regions {wildcards.part} --output-type b {input.vcf} > {output}
        """

#TODO investigate filter option
rule STITCH_format:
    input: "Imputation/input/{part}.bisnp.bcf"
    output: "Imputation/input/{part}.stitch"
    message: "Converting biallelic data to STITCH format: {wildcards.part}"
    benchmark: "Benchmark/Impute/stitchformat.{part}.txt"
    threads: 1
    params: 
        filters = "-i \'QUAL>20 && DP>10\'" if config["filtervcf"] else ""
    shell:
        """
        bcftools query {params} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input} > {output}
        """

rule impute:
    input:
        bamlist = "Imputation/input/samples.list",
        infile = "Imputation/input/{part}.stitch",
        chromosome = "Imputation/input/contigs/{part}"
    output:
        # format a wildcard pattern like "k{k}/s{s}/ngen{ngen}"
        # into a file path, with k, s, ngen being the columns of the data frame
        f"Imputation/{paramspace.wildcard_pattern}/contigs/" + "{part}/impute.vcf.gz"
    log: f"Imputation/{paramspace.wildcard_pattern}/contigs/" + "{part}/stitch.log"
    params:
        # automatically translate the wildcard values into an instance of the param space
        # in the form of a dict (here: {"k": ..., "s": ..., "ngen": ...})
        parameters = paramspace.instance
    message: "Running STITCH: {wildcards.part}\n  Parameters:\n  " + "{params.parameters}"
    benchmark: f"Benchmark/Impute/stitch.{paramspace.wildcard_pattern}" + ".{part}.txt"
    threads: 50
    script: "../utilities/stitch_impute.R"

rule samples_file:
    output: "Imputation/input/samples.names"
    message: "Creating file of sample names"
    threads: 1
    run:
        with open(output[0], "w") as fout:
            [fout.write(f"{i}\n") for i in samplenames]

rule index_vcf:
    input:
        vcf = "Imputation/{stitchparams}/contigs/{part}/impute.vcf.gz",
        samplelist = "Imputation/input/samples.names"
    output: "Imputation/{stitchparams}/contigs/{part}/impute.vcf.gz.tbi"
    log: "Imputation/{stitchparams}/contigs/{part}/impute.stats"
    message: "Indexing: {wildcards.stitchparams}/{wildcards.part}"
    benchmark: "Benchmark/Impute/indexvcf.{stitchparams}.{part}.txt"
    threads: 1
    shell:
        """
        tabix {input.vcf}
        bcftools stats {input.vcf} -S {input.samplelist} > {log}
        """

rule merge_vcfs:
    input: 
        vcf = expand("Imputation/{{stitchparams}}/contigs/{part}/impute.vcf.gz", part = contigs),
        idx = expand("Imputation/{{stitchparams}}/contigs/{part}/impute.vcf.gz.tbi", part = contigs),
    output: "Imputation/{stitchparams}/variants.imputed.bcf"
    log: "Imputation/{stitchparams}/concat.log"
    message: "Merging VCFs: {wildcards.stitchparams}"
    benchmark: "Benchmark/Impute/mergevcf.{stitchparams}.txt"
    threads: 20
    shell:
        """
        bcftools concat --threads {threads} -o {output} --output-type b {input.vcf} 2> {log}
        """

rule stats:
    input:
        bcf = "Imputation/{stitchparams}/variants.imputed.bcf",
        samplelist = "Imputation/input/samples.names"
    output: "Imputation/{stitchparams}/variants.imputed.stats"
    message: "Indexing and calculating stats: {wildcards.stitchparams}/variants.imputed.bcf"
    benchmark: "Benchmark/Impute/mergestats.{stitchparams}.txt"
    shell:
        """
        bcftools index {input.bcf}
        bcftools stats {input.bcf} -S {input.samplelist} > {output}
        """

rule reports:
    input: "Imputation/{stitchparams}/variants.imputed.stats"
    output: "Imputation/{stitchparams}/variants.imputed.html"
    message: "Generating bcftools report: {output}"
    benchmark: "Benchmark/Impute/stitchreport.{stitchparams}.txt"
    script: "../utilities/bcftoolsreport.Rmd"

rule all:
    input: 
        bcf = expand("Imputation/{stitchparams}/variants.imputed.bcf", stitchparams=paramspace.instance_patterns),
        reports = expand("Imputation/{stitchparams}/variants.imputed.html", stitchparams=paramspace.instance_patterns)
    default_target: True