import os

bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
groupings = config.get("groupings", None)
ploidy = config["ploidy"]
samplenames = config["samplenames"]
mp_extra = config.get("extra", "") 

if groupings is not None:
	absent = []
	with open(groupings) as f:
		for line in f:
			samp, pop = line.rstrip().split()
			if samp not in samplenames:
				absent.append(samp)
	if absent:
		sys.tracebacklimit = 0
		raise ValueError(f"{len(absent)} sample(s) in \033[1m{groupings}\033[0m not found in \033[1m{bam_dir}\033[0m directory:\n\033[33m" + ", ".join(absent) + "\033[0m")

def faidx_contignames(infile):
	bn = os.path.basename(infile)
	os.makedirs("Assembly", exist_ok = True)
	if not os.path.exists(f"Assembly/{bn}"):
		shell(f"ln -sr {infile} Assembly/{bn}")
	if not os.path.exists(f"Assembly/{bn}.fai"):
		print(f"Assembly/{bn}.fai not found, indexing {bn} with samtools faidx")
		subprocess.run(["samtools","faidx", "--fai-idx", f"Assembly/{bn}.fai", infile, "2>", "/dev/null"])
	with open(f"Assembly/{bn}.fai") as f:
		lines = [line.rstrip().split("\t")[0] for line in f]
	return lines

contigs = faidx_contignames(genomefile)
dict_cont = dict(zip(contigs, contigs))

rule index_alignments:
	input: bam_dir + "/{sample}.bam"
	output: bam_dir + "/{sample}.bam.bai"
	message: "Indexing barcodes: {wildcards.sample}"
	benchmark: "Benchmark/Variants/mpileup/indexbam.{sample}.txt"
	shell:
		"sambamba index {input} {output}"

rule bam_list:
	input: 
		bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
		bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames)
	output: "Variants/mpileup/logs/samples.files"
	message: "Creating list of alignment files"
	benchmark: "Benchmark/Variants/mpileup/bamlist.txt"
	run:
		with open(output[0], "w") as fout:
			for bamfile in input.bam:
				_ = fout.write(bamfile + "\n")

rule samplenames:
	output: "Variants/mpileup/logs/samples.names"
	message: "Creating list of sample names"
	run:
		with open(output[0], "w") as fout:
			for samplename in samplenames:
				_ = fout.write(samplename + "\n")		

rule mpileup:
	input:
		bamlist = "Variants/mpileup/logs/samples.files",
		genome = f"Assembly/{genomefile}"
	output: 
		pipe("Variants/mpileup/{part}.mp.bcf")
	params: 
		lambda wc: dict_cont[wc.part]
	message: 
		"Finding variants: {wildcards.part}"
	log: 
		"Variants/mpileup/logs/{part}.mpileup.log"
	benchmark: 
		"Benchmark/Variants/mpileup/mpileup.{part}.txt"
	params:
		region = "{wildcards.part}",
		extra = mp_extra
	shell:
		"bcftools mpileup --fasta-ref {input.genome} --region {params} --bam-list {input.bamlist} --annotate AD --output-type b > {output} 2> {log}"

rule call_genotypes:
	input: "Variants/mpileup/{part}.mp.bcf"
	output: temp("Variants/mpileup/call/{part}.bcf")
	message: "Calling genotypes: {wildcards.part}"
	benchmark: "Benchmark/Variants/mpileup/call.{part}.txt"
	log: "Variants/mpileup/logs/{part}.call.log"
	threads: 2
	params: 
		groupsamples = '' if groupings is None else f"--group-samples {groupings}",
		ploidy = f"--ploidy {ploidy}"
	shell:
		"bcftools call --multiallelic-caller {params} --variants-only --output-type b {input} | bcftools sort - --output {output} 2> /dev/null"

rule index_bcf:
	input: 
		bcf = "Variants/mpileup/call/{part}.bcf",
		samplelist = "Variants/mpileup/logs/samples.names",
		genome = f"Assembly/{genomefile}"
	output: temp("Variants/mpileup/call/{part}.bcf.csi")
	log: "Variants/mpileup/stats/{part}.stats"
	message: "Indexing: {wildcards.part}"
	benchmark: "Benchmark/Variants/mpileup/indexbcf.{part}.txt"
	threads: 4
	shell:
		"""
		bcftools index --threads {threads} --output {output} {input.bcf}
		bcftools stats -S {input.samplelist} --fasta-ref {input.genome} {input.bcf} > {log}
		"""

rule combine_bcfs:
	input: 
		bcf = expand("Variants/mpileup/call/{part}.bcf", part = contigs),
		idx = expand("Variants/mpileup/call/{part}.bcf.csi", part = contigs),
		genome = f"Assembly/{genomefile}",
		samplelist = "Variants/mpileup/logs/samples.names"
	output: 
		bcf = "Variants/mpileup/variants.raw.bcf",
		idx = "Variants/mpileup/variants.raw.bcf.csi",
		stats = "Variants/mpileup/stats/variants.raw.stats"
	message: "Merging all BCFs into: {output.bcf}"
	benchmark: "Benchmark/Variants/mpileup/merge.txt"
	threads: 50
	shell:
		"""
		bcftools concat --threads {threads} --output-type b --naive {input.bcf} > {output.bcf} 2> /dev/null
		bcftools index --output {output.idx} {output.bcf}
		bcftools stats -S {input.samplelist} --fasta-ref {input.genome} {output.bcf} > {output.stats}
		"""

rule normalize_bcf:
	input: 
		genome = f"Assembly/{genomefile}",
		bcf = "Variants/mpileup/variants.raw.bcf",
		samplelist = "Variants/mpileup/logs/samples.names"
	output:
		bcf = "Variants/mpileup/variants.normalized.bcf",
		idx = "Variants/mpileup/variants.normalized.bcf.csi",
		stats = "Variants/mpileup/stats/variants.normalized.stats"
	message: 
		"Normalizing the called variants"
	threads: 2
	shell:
		"""
		bcftools norm -d none -f {input.genome} {input.bcf} | bcftools norm -m -any -N -Ob > {output.bcf}
		bcftools index --output {output.idx} {output.bcf}
		bcftools stats -S {input.samplelist} --fasta-ref {input.genome} {output.bcf} > {output.stats}
		"""

rule bcfreport:
	input: "Variants/mpileup/stats/variants.raw.stats"
	output: "Variants/mpileup/stats/variants.raw.html"
	message: "Generating bcftools report: variants.raw.bcf"
	benchmark: "Benchmark/Variants/mpileup/reports.txt"
	script: "../utilities/reportBcftools.Rmd"

rule bcfreportnorm:
	input: "Variants/mpileup/stats/variants.normalized.stats"
	output: "Variants/mpileup/stats/variants.normalized.html"
	message: "Generating bcftools report: variants.normalized.bcf"
	script: "../utilities/reportBcftools.Rmd"

rule all:
	input: 
		expand("Variants/mpileup/variants.{file}.bcf", file = ["raw","normalized"]),
		expand("Variants/mpileup/stats/variants.{file}.html", file = ["raw","normalized"])
	default_target: True
	message: "Variant calling is complete!"