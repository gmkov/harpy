import os

bam_dir     = config["seq_directory"]
genomefile  = config["genomefile"]
samplenames = config["samplenames"] 
extra       = config.get("extra", "") 
outdir      = "Variants/leviathan"
bn          = os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
if genome_zip:
    bn = bn[:-3]

conda:
    os.getcwd() + "/harpyenvs/variants.sv.yaml"

rule index_alignment:
    input:
        bam_dir + "/{sample}.bam"
    output:
        bam_dir + "/{sample}.bam.bai"
    message:
        "Indexing alignment: {wildcards.sample}"
    shell:
        "samtools index {input} {output} 2> /dev/null"

rule index_barcode:
    input: 
        bam = bam_dir + "/{sample}.bam",
        bai = bam_dir + "/{sample}.bam.bai"
    output:
        temp(outdir + "/lrezIndexed/{sample}.bci")
    message:
        "Indexing barcodes: {wildcards.sample}"
    benchmark:
        ".Benchmark/Variants/leviathan/indexbc.{sample}.txt"
    threads:
        4
    shell:
        "LRez index bam --threads {threads} -p -b {input.bam} -o {output}"

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    message: 
        "Creating {output}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be decompressed
            gzip -dc {input} > {output}
        elif (file {input} | grep -q BGZF ); then
            # is bgzipped, decompress
            gzip -dc {input} > {output}
        else
            # isn't compressed, just linked
            ln -sr {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai"
    message:
        "Indexing {input}"
    log:
        f"Genome/{bn}.faidx.log"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule index_bwa_genome:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    message:
        "Indexing {input}"
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    log:
        f"Genome/{bn}.idx.log"
    shell: 
        "bwa index {input} 2> {log}"

rule leviathan_variantcall:
    input:
        bam    = bam_dir + "/{sample}.bam",
        bai    = bam_dir + "/{sample}.bam.bai",
        bc_idx = outdir + "/lrezIndexed/{sample}.bci",
        genome = f"Genome/{bn}",
        genidx = multiext(f"Genome/{bn}", ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        pipe(outdir + "/{sample}.vcf")
    log:  
        runlog     = outdir + "/logs/{sample}.leviathan.log",
        candidates = outdir + "/logs/{sample}.candidates"
    message:
        "Calling variants: {wildcards.sample}"
    benchmark:
        ".Benchmark/Variants/leviathan/variantcall.{sample}.txt"
    params:
        extra = extra
    threads:
        3
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}"

rule sort_bcf:
    input:
        outdir + "/{sample}.vcf"
    output:
        outdir + "/{sample}.bcf"
    message:
        "Sorting and converting to BCF: {wildcards.sample}"
    params:
        "{wildcards.sample}"
    benchmark:
        ".Benchmark/Variants/leviathan/sortbcf.{sample}.txt"
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule sv_stats:
    input: 
        outdir + "/{sample}.bcf"
    output: 
        outdir + "/reports/stats/{sample}.sv.stats"
    message:
        "Getting SV stats for {wildcards.sample}"
    benchmark:
        ".Benchmark/Variants/leviathan/stats.{sample}.txt"
    shell:
        """
        echo -e "sample\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.sample}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule sv_report:
    input:	
        bcf       = outdir + "/{sample}.bcf",
        statsfile = outdir + "/reports/stats/{sample}.sv.stats"
    output:	
        outdir + "/reports/{sample}.SV.html"
    message:
        "Generating SV report: {wildcards.sample}"
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    script:
        "reportLeviathan.Rmd"

rule log_runtime:
    output:
        outdir + "/logs/harpy.variants.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        extra = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy variants sv module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write("The barcodes were indexed using:\n")
            _ = f.write("    LRez index bam -p -b INPUT\n")
            _ = f.write("Leviathan was called using:\n")
            _ = f.write(f"    LEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}\n")

rule all_bcfs:
    default_target: True
    input: 
        bcf     = expand(outdir + "/{sample}.bcf", sample = samplenames),
        reports = expand(outdir + "/reports/{sample}.SV.html", sample = samplenames),
        runlog  = outdir + "/logs/harpy.variants.log"
    message:
        "Variant calling is complete!"
