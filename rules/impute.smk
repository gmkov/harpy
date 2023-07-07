from snakemake.utils import Paramspace
import pandas as pd
import subprocess
import sys

bam_dir     = config["seq_directory"]
samplenames = config["samplenames"]
variantfile = config["variantfile"]
paramfile   = config["paramfile"]
paramspace  = Paramspace(pd.read_csv(paramfile, sep="\t"), param_sep = "", filename_params="*")
#^ declare a dataframe to be a paramspace

def contignames(vcf):
    print("Preprocessing: Identifying contigs with at least 2 biallelic SNPs", file = sys.stderr)
    biallelic = subprocess.Popen(f"bcftools view -m2 -M2 -v snps {vcf} -Ob".split(), stdout = subprocess.PIPE)
    contigs = subprocess.run("""bcftools query -f %CHROM\\n""".split(), stdin = biallelic.stdout, stdout = subprocess.PIPE)
    dict_cont = dict()
    for i in list([chr for chr in contigs.stdout.decode('utf-8').split()]):
        if i in dict_cont:
            dict_cont[i] += 1
        else:
            dict_cont[i] = 1
    return [contig for contig in dict_cont if dict_cont[contig] > 1]

contigs   = contignames(variantfile)
dict_cont = dict(zip(contigs, contigs))

rule bam_list:
    input:
        expand(bam_dir + "/{sample}.bam", sample = samplenames)
    output:
        "Impute/input/samples.list"
    message:
        "Creating list of alignment files"
    benchmark:
        "Benchmark/Impute/filelist.txt"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input:
                fout.write(f"{bamfile}\n")

rule samples_file:
    output:
        "Impute/input/samples.names"
    message:
        "Creating file of sample names"
    threads: 1
    run:
        with open(output[0], "w") as fout:
            [fout.write(f"{i}\n") for i in samplenames]

##TODO investigate filter option

rule sort_bcf:
    input:
        variantfile
    output:
        bcf = temp("Impute/input/input.sorted.bcf"),
        idx = temp("Impute/input/input.sorted.bcf.csi")
    log:
        "Impute/input.sorted.log"
    message:
        "Sorting input variant call file"
    shell:
        """
        bcftools sort -Ob {input} > {output.bcf} 2> {log}
		bcftools index --output {output.idx} {output.bcf}
        """

rule convert2stitch:
    input:
        variantfile
    output:
        "Impute/input/{part}.stitch"
    message:
        "Converting data to biallelic STITCH format: {wildcards.part}"
    #params:
        #filters = "-i \'QUAL>20 && DP>10\'" if config["filtervcf"] else ""
    benchmark:
        "Benchmark/Impute/fileprep.{part}.txt"
    threads: 3
    shell:
        """
        bcftools view -m2 -M2 -v snps --regions {wildcards.part} {input} |
        bcftools query -i 'TYPE="SNP"' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' > {output}
        """

rule impute:
    input:
        bamlist = "Impute/input/samples.list",
        infile  = "Impute/input/{part}.stitch"
    output:
        # format a wildcard pattern like "k{k}/s{s}/ngen{ngen}"
        # into a file path, with k, s, ngen being the columns of the data frame
        f"Impute/{paramspace.wildcard_pattern}/contigs/" + "{part}/{part}.vcf.gz"
    log:
        f"Impute/{paramspace.wildcard_pattern}/contigs/" + "{part}/{part}.log"
    params:
        # automatically translate the wildcard values into an instance of the param space
        # in the form of a dict (here: {"k": ..., "s": ..., "ngen": ...})
        parameters = paramspace.instance
    message: 
        "Running STITCH: {wildcards.part}\nmodel: {wildcards.model}\nuseBX: {wildcards.useBX}\n    k: {wildcards.k}\n    s: {wildcards.s}\n nGen: {wildcards.nGen}"
    benchmark:
        f"Benchmark/Impute/stitch.{paramspace.wildcard_pattern}" + ".{part}.txt"
    threads: 50
    script:
        "stitch_impute.R"


rule index_vcf:
    input:
        vcf        = "Impute/{stitchparams}/contigs/{part}/{part}.vcf.gz",
        samplelist = "Impute/input/samples.names"
    output: 
        idx        = "Impute/{stitchparams}/contigs/{part}/{part}.vcf.gz.tbi",
        stats      = "Impute/{stitchparams}/contigs/{part}/{part}.stats"
    message:
        "Indexing: {wildcards.stitchparams}/{wildcards.part}"
    benchmark:
        "Benchmark/Impute/indexvcf.{stitchparams}.{part}.txt"
    threads: 1
    shell:
        """
        tabix {input.vcf}
        bcftools stats {input.vcf} -S {input.samplelist} > {output.stats}
        """

rule stitch_reports:
    input:
        "Impute/{stitchparams}/contigs/{part}/{part}.stats"
    output:
        "Impute/{stitchparams}/contigs/{part}/{part}.impute.html"
    message:
        "Generating STITCH report: {wildcards.part}"
    benchmark:
        "Benchmark/Impute/report.{stitchparams}.{part}.txt"
    threads: 1
    script:
        "reportStitch.Rmd"

rule clean_stitch:
    input:
        "Impute/{stitchparams}/contigs/{part}/{part}.impute.html"
    output:
        temp("Impute/{stitchparams}/contigs/{part}/.cleaned")
    message:
        "Cleaning up {wildcards.stitchparams}: {wildcards.part}"
    priority: 1
    shell: 
        """
        rm -rf Impute/{wildcards.stitchparams}/contigs/{wildcards.part}/input
        rm -rf Impute/{wildcards.stitchparams}/contigs/{wildcards.part}/RData
        rm -rf Impute/{wildcards.stitchparams}/contigs/{wildcards.part}/plots
        touch {output}
        """

rule merge_vcfs:
    input: 
        vcf   = expand("Impute/{{stitchparams}}/contigs/{part}/{part}.vcf.gz", part = contigs),
        idx   = expand("Impute/{{stitchparams}}/contigs/{part}/{part}.vcf.gz.tbi", part = contigs),
        clean = expand("Impute/{{stitchparams}}/contigs/{part}/.cleaned", part = contigs)
    output:
        "Impute/{stitchparams}/variants.imputed.bcf"
    log:
        "Impute/{stitchparams}/concat.log"
    message:
        "Merging VCFs: {wildcards.stitchparams}"
    benchmark:
        "Benchmark/Impute/mergevcf.{stitchparams}.txt"
    threads: 50
    shell:
        "bcftools concat --threads {threads} -o {output} --output-type b {input.vcf} 2> {log}"

rule stats:
    input:
        bcf        = "Impute/{stitchparams}/variants.imputed.bcf",
        samplelist = "Impute/input/samples.names"
    output:
        "Impute/{stitchparams}/variants.imputed.stats"
    message:
        "Indexing and calculating stats: {wildcards.stitchparams}/variants.imputed.bcf"
    benchmark:
        "Benchmark/Impute/mergestats.{stitchparams}.txt"
    shell:
        """
        bcftools index {input.bcf}
        bcftools stats {input.bcf} -S {input.samplelist} > {output}
        """

rule reports:
    input: 
        "Impute/{stitchparams}/variants.imputed.stats"
    output:
        "Impute/{stitchparams}/variants.imputed.html"
    message:
        "Generating bcftools report: {output}"
    benchmark:
        "Benchmark/Impute/stitchreport.{stitchparams}.txt"
    script:
        "reportBcftools.Rmd"

rule log_runtime:
    output:
        "Impute/logs/harpy.impute.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy impute module ran using these parameters:\n\n")
            _ = f.write(f"The provided variant file: {variantfile}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write("Preprocessing was performed with:\n")
            _ = f.write("\tbcftools view -m2 -M2 -v snps --regions CONTIG INFILE |\n")
            _ = f.write("\tbcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'\n")
            _ = f.write("## STITCH imputation ##\n")
            _ = f.write("The STITCH parameters were governed by the rows of the input parameter table:\n")
            with open(config["paramfile"], "r") as f1:
                for line in f1:
                    _ = f.write(line)
            _ = f.write(
                """
                STITCH(
                    method               = model,
                    posfile              = posfile,
                    bamlist              = bamlist,
                    nCores               = nCores,
                    nGen                 = nGen,
                    chr                  = chr,
                    K                    = k,
                    S                    = s,
                    use_bx_tag           = useBX,
                    bxTagUpperLimit      = 50000,
                    niterations          = 40,
                    switchModelIteration = 39,
                    splitReadIterations  = NA,
                    outputdir            = outdir,
                    output_filename      = outfile
                )
                """
            )

rule all:
    input: 
        bcf           = expand("Impute/{stitchparams}/variants.imputed.bcf", stitchparams=paramspace.instance_patterns),
        reports       = expand("Impute/{stitchparams}/variants.imputed.html", stitchparams=paramspace.instance_patterns),
        contigreports = expand("Impute/{stitchparams}/contigs/{part}/{part}.impute.html", stitchparams=paramspace.instance_patterns, part = contigs),
        runlog        = "Impute/logs/harpy.impute.log"
    message: 
        "Genotype imputation is complete!"
    default_target: True