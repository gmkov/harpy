import os
import re

seq_dir 	= config["seq_directory"]
nbins 		= config["EMA_bins"]
genomefile 	= config["genomefile"]
Rsep 		= config["Rsep"]
fqext 		= config["fqext"]
#samplenames = config["samplenames"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
outdir      = "Align/ema"

flist = os.listdir(seq_dir)
r = re.compile(".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
fqlist = list(filter(r.match, flist))
bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

rule all:
    input: 
        expand(outdir + "/{sample}.bam", sample = samplenames),
        expand(outdir + "/{sample}.bam.bai", sample = samplenames),
        expand(outdir + "/stats/BXstats/{sample}.bxstats.html", sample = samplenames),
        expand(outdir + "/stats/coverage/{sample}.cov.html", sample = samplenames),
        outdir + "/stats/reads.bxcounts.html",
        outdir + "/stats/samtools_stats/alignment.stats.html",
        outdir + "/stats/samtools_flagstat/alignment.flagstat.html",
        outdir + "/logs/harpy.align.log"
    message:
        "Read mapping completed!"
    benchmark:
        "Benchmark/Mapping/ema/report.txt"
    default_target: True

rule link_genome:
    input:
        genomefile
    output: 
        f"Assembly/{bn}"
    message:
        "Symlinking {input} to Assembly/"
    shell: 
        "ln -sr {input} {output}"


rule faidx_genome:
    input: 
        f"Assembly/{bn}"
    output: 
        f"Assembly/{bn}.fai"
    message:
        "Indexing {input}"
    log:
        f"Assembly/{bn}.faidx.log"
    shell: 
        """
        samtools faidx --fai-idx {output} {input} 2> {log}
        """

rule index_bwa_genome:
    input: 
        f"Assembly/{bn}"
    output: 
        multiext(f"Assembly/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    message:
        "Indexing {input}"
    log:
        f"Assembly/{bn}.idx.log"
    shell: 
        """
        bwa index {input} 2> {log}
        """

rule make_genome_windows:
    input:
        f"Assembly/{bn}.fai"
    output: 
        f"Assembly/{bn}.bed"
    message: 
        "Creating BED intervals from {input}"
    shell: 
        """
        makewindows.py -i {input} -w 10000 -o {output}
        """

rule count_beadtags:
    input:
        forward_reads = seq_dir + "/{sample}" + f"{Rsep[0]}.{fqext}",
        reverse_reads = seq_dir + "/{sample}" + f"{Rsep[1]}.{fqext}"
    output: 
        counts = outdir + "/count/{sample}.ema-ncnt",
        logs   = temp(outdir + "/count/logs/{sample}.count.log")
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    message:
        "Counting barcode frequency: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/ema/Count.{sample}.txt"
    params:
        prefix = lambda wc: outdir + "/count/" + wc.get("sample")
    threads: 1
    shell:
        "seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema count -p -o {params} 2> {output.logs}"

rule beadtag_summary:
    input: 
        countlog = expand(outdir + "/count/logs/{sample}.count.log", sample = samplenames)
    output:
        outdir + "/stats/reads.bxcounts.html"
    message:
        "Creating sample barcode validation report"
    benchmark:
        "Benchmark/Mapping/ema/beadtagsummary.txt"
    script:
        "reportEmaCount.Rmd"

rule preprocess_ema:
    input: 
        forward_reads = seq_dir + "/{sample}" + f"{Rsep[0]}.{fqext}",
        reverse_reads = seq_dir + "/{sample}" + f"{Rsep[1]}.{fqext}",
        emacounts     = outdir + "/count/{sample}.ema-ncnt"
    output: 
        bins       	  = temp(expand(outdir + "/preproc/{{sample}}/ema-bin-{bin}", bin = ["%03d" % i for i in range(nbins)])),
        unbarcoded    = temp(outdir + "/preproc/{sample}/ema-nobc")
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    log:
        outdir + "/preproc/logs/{sample}.preproc.log"
    message:
        "Preprocessing for EMA mapping: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/ema/Preproc.{sample}.txt"
    threads:
        2
    params:
        outdir = lambda wc: outdir + "/preproc/" + wc.get("sample"),
        bins   = nbins
    shell:
        "seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema preproc -p -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 | cat - > {log}"

rule align_ema:
    input:
        readbin    = outdir + "/preproc/{sample}/ema-bin-{bin}",
        genome 	   = f"Assembly/{bn}",
        genome_idx = multiext(f"Assembly/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    output:
        alignment  = temp(outdir + "/align/{sample}/{sample}.{bin}.bam")
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    message:
        "Aligning barcoded sequences: {wildcards.sample}-{wildcards.bin}"
    benchmark:
        "Benchmark/Mapping/ema/Align.{sample}.{bin}.txt"
    params: 
        quality = config["quality"],
        extra = extra
    threads: 8
    shell:
        """
        EMATHREADS=$(( {threads} - 2 ))
        ema align -t $EMATHREADS {params.extra} -d -p haptag -r {input.genome} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" -s {input.readbin} 2> /dev/null |
        samtools view -h -F 4 -q {params.quality} - | 
        samtools sort --reference {input.genome} -O bam -m 4G -o {output} - 2> /dev/null
        """

rule align_nobarcode:
    input:
        reads      = outdir + "/preproc/{sample}/ema-nobc",
        genome 	   = f"Assembly/{bn}",
        genome_idx = multiext(f"Assembly/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    output: 
        samfile    = temp(outdir + "/align/{sample}/{sample}.nobarcode.bam.tmp")
    benchmark:
        "Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    params:
        quality = config["quality"]
    message:
        "Aligning unbarcoded sequences: {wildcards.sample}"
    threads: 8
    shell:
        """
        BWATHREADS=$(( {threads} - 2 ))
        bwa mem -t $BWATHREADS -C -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.reads} 2> /dev/null |
        samtools view -h -F 4 -q {params.quality} | 
        samtools sort -O bam -m 4G --reference {input.genome} -o {output} 2> /dev/null
        """

rule markduplicates:
    input:
        bam      = outdir + "/align/{sample}/{sample}.nobarcode.bam.tmp"
    output: 
        bam      = temp(outdir + "/align/{sample}/{sample}.nobarcode.bam"),
        bai      = temp(outdir + "/align/{sample}/{sample}.nobarcode.bam.bai")
    log: 
        mdlog    = outdir + "/stats/markduplicates/{sample}.markdup.nobarcode.log",
        stats    = outdir + "/stats/samtools_stats/{sample}.nobarcode.stats",
        flagstat = outdir + "/stats/samtools_flagstat/{sample}.nobarcode.flagstat"
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    message:
        "Marking duplicates in unbarcoded alignments: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/ema/markdup.{sample}.txt"
    threads: 2
    shell:
        """
        sambamba markdup -t {threads} -l 4 {input} {output.bam} 2> {log.mdlog}
        samtools stats {output.bam} > {log.stats}
        samtools flagstat {output.bam} > {log.flagstat}
        """   

rule merge_barcoded:
    input:
        aln = expand(outdir + "/align/{{sample}}/{{sample}}.{bin}.bam", bin = ["%03d" % i for i in range(nbins)]),
    output: 
        bam = temp(outdir + "/align/barcoded/{sample}.barcoded.bam"),
        bai = temp(outdir + "/align/barcoded/{sample}.barcoded.bam.bai")
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    message:
        "Merging barcoded alignments: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/ema/merge.{sample}.txt"
    threads: 10
    shell:
        "sambamba merge -t {threads} -l 4 {output.bam} {input} 2> /dev/null"

#rule secondary2split:
#	input:
#		bam = outdir + "/align/barcoded/{sample}.barcoded.sec.bam",
#		bai = outdir + "/align/barcoded/{sample}.barcoded.sec.bam.bai"
#	output:
#		bam = temp(outdir + "/align/barcoded/{sample}.barcoded.bam"),
#		bai = temp(outdir + "/align/barcoded/{sample}.barcoded.bam.bai")
#	wildcard_constraints:
#		sample = "[a-zA-Z0-9\_\-\.]*"
#	message:
#		"Converting Secondary SAM flags to Split flags: {wildcards.sample}"
#	shell:
#		"secondary2split.py {input.bam} {output.bam}"

rule bcstats:
    input: 
        bam      = outdir + "/align/barcoded/{sample}.barcoded.bam",
        bai      = outdir + "/align/barcoded/{sample}.barcoded.bam.bai"
    log:
        stats    = outdir + "/stats/samtools_stats/{sample}.barcoded.stats",
        flagstat = outdir + "/stats/samtools_flagstat/{sample}.barcoded.flagstat"
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    message:
        "Indexing merged barcoded alignemnts: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/ema/indexmerge.{sample}.txt"
    shell:
        """
        samtools stats {input.bam} > {log.stats}
        samtools flagstat {input.bam} > {log.flagstat}
        """

rule alignment_coverage:
    input: 
        bed     = f"Assembly/{bn}.bed",
        nobx    = outdir + "/align/{sample}/{sample}.nobarcode.bam",
        nobxbai = outdir + "/align/{sample}/{sample}.nobarcode.bam.bai",
        bx      = outdir + "/align/barcoded/{sample}.barcoded.bam",
        bxbai   = outdir + "/align/barcoded/{sample}.barcoded.bam.bai"
    output: 
        outdir + "/stats/coverage/data/{sample}.cov.gz"
    message:
        "Calculating genomic coverage: {wildcards.sample}"
    threads: 2
    shell:
        "samtools bedcov -c {input.bed} {input.bx} {input.nobx} | gzip > {output}"

rule gencovBX_report:
    input: 
        outdir + "/stats/coverage/data/{sample}.cov.gz",
    output:
        outdir + "/stats/coverage/{sample}.cov.html"
    message:
        "Creating report of alignment coverage: {wildcards.sample}"
    script:
        "reportEmaGencov.Rmd"

rule merge_alignments:
    input:
        aln_barcoded  = outdir + "/align/barcoded/{sample}.barcoded.bam",
        idx_barcoded  = outdir + "/align/barcoded/{sample}.barcoded.bam.bai",
        aln_nobarcode = outdir + "/align/{sample}/{sample}.nobarcode.bam",
        idx_nobarcode = outdir + "/align/{sample}/{sample}.nobarcode.bam.bai"
    output: 
        bam 		  = temp(outdir + "/align/{sample}.unsort.bam"),
        bai 		  = temp(outdir + "/align/{sample}.unsort.bam.bai")
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    message:
        "Merging all alignments: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/ema/mergebc_nobc.{sample}.txt"
    threads: 10
    shell:
        "sambamba merge -t {threads} {output.bam} {input.aln_barcoded} {input.aln_nobarcode} 2> /dev/null"

rule sort_merge:
    input:
        bam    = outdir + "/align/{sample}.unsort.bam",
        genome = f"Assembly/{bn}"
    output:
        outdir + "/{sample}.bam"
    message:
        "Sorting merged barcoded alignments: {wildcards.sample}"
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    threads: 2
    priority: 1
    shell:
        "samtools sort -@ {threads} -O bam --reference {input.genome} -m 4G -o {output} {input.bam} 2> /dev/null"

rule index_alignments:
    input: 
        outdir + "/{sample}.bam"
    output:
        outdir + "/{sample}.bam.bai"
    message:
        "Indexing: {input}"
    benchmark:
        "Benchmark/Mapping/ema/IndexMerged.{sample}.txt"
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    shell:
        "sambamba index {input} {output} 2> /dev/null"

rule alignment_bxstats:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/stats/BXstats/data/{sample}.bxstats.gz"
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    threads: 1
    shell:
        "bxStats.py {input.bam} > {output}"

rule bx_stats_report:
    input:
        outdir + "/stats/BXstats/data/{sample}.bxstats.gz"
    output:	
        outdir + "/stats/BXstats/{sample}.bxstats.html"
    message: 
        "Generating summary of barcode alignment: {wildcards.sample}"
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    threads: 1
    script:
        "reportBxStats.Rmd"

rule general_alignment_stats:
    input: 		
        bam      = outdir + "/{sample}.bam",
        bai      = outdir + "/{sample}.bam.bai"
    output:
        stats    = outdir + "/stats/samtools_stats/{sample}.stats",
        flagstat = outdir + "/stats/samtools_flagstat/{sample}.flagstat"
    message:
        "Calculating alignment stats: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/ema/Mergedstats.{sample}.txt"
    wildcard_constraints:
        sample = "[a-zA-Z0-9\_\-\.]*"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_reports:
    input: 
        expand(outdir + "/stats/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
    output: 
        stats    = outdir + "/stats/samtools_stats/alignment.stats.html",
        flagstat = outdir + "/stats/samtools_flagstat/alignment.flagstat.html"
    message:
        "Summarizing samtools stats and flagstats"
    benchmark:
        "Benchmark/Mapping/ema/report.txt"
    shell:
        """
        multiqc Align/ema/stats/samtools_stats    --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
        multiqc Align/ema/stats/samtools_flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
        """

rule log_runtime:
    output:
        outdir + "/logs/harpy.align.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        extra = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy align module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("Barcodes were counted and validated with EMA using:\n")
            _ = f.write("\tseqfu interleave forward.fq.gz reverse.fq.gz | ema count -p\n")
            _ = f.write("Barcoded sequences were binned with EMA using:\n")
            _ = f.write(f"\tseqfu interleave forward.fq.gz reverse.fq.gz | ema preproc -p -n {nbins}\n")
            _ = f.write("Barcoded bins were aligned with ema align using:\n")
            _ = f.write("\tema align " + extra + " -d -p haptag -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" |\n")
            _ = f.write("\tsamtools view -h -F 4 -q " + str(config["quality"]) + " - |\n") 
            _ = f.write("\tsamtools sort --reference genome -m 4G\n\n")
            _ = f.write("Invalid/non barcoded sequences were aligned with BWA using:\n")
            _ = f.write("\tbwa mem -C -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome forward_reads reverse_reads |\n")
            _ = f.write("\tsambamba markdup -l 0\n")
