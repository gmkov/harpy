from .helperfunctions import generate_conda_deps, fetch_file 
from .printfunctions import print_onstart
from .validations import validate_input_by_ext
import rich_click as click
import subprocess
import os
import sys

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/simulate")
@click.option('-d', '--outer-distance', type = click.IntRange(min = 100), default = 350, show_default= True, help = "Outer distance between paired-end reads (bp)")
@click.option('-i', '--distance-sd', type = click.IntRange(min = 1), default = 15, show_default=True,  help = "Standard deviation of read-pair distance")
@click.option('-b', '--barcodes', type = click.Path(exists=True, dir_okay=False), help = "File of linked-read barcodes")
@click.option('-n', '--read-pairs', type = click.IntRange(min = 1), default = 600, show_default=True,  help = "Number of read pairs to simulate, in millions")
@click.option('-l', '--molecule-length', type = click.IntRange(min = 10), default = 100, show_default=True,  help = "Mean molecule length (kbp)")
@click.option('-p', '--partitions', type = click.IntRange(min = 1), default=1500, show_default=True,  help = "How many partitions to generate (×1000)")
@click.option('-m', '--molecules-per', type = click.IntRange(min = 1), default = 100, show_default=True,  help = "Average number of molecules per partition")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Print the generated snakemake command and exit')
@click.option('-o', '--output-dir', type = str, default = "Simulate/linkedreads", help = 'Name of output directory')
@click.argument('genome_hap1', required=True, type=click.Path(exists=True, dir_okay=False), nargs=1)
@click.argument('genome_hap2', required=True, type=click.Path(exists=True, dir_okay=False), nargs=1)
def linkedreads(genome_hap1, genome_hap2, output_dir, outer_distance, distance_sd, barcodes, read_pairs, molecule_length, partitions, molecules_per, threads, snakemake, quiet, print_only):
    """
    Create linked reads from a genome
 
    If not providing a text file of `--barcodes` to use for the simulated linked reads, Harpy will
    download the `4M-with-alts-february-2016.txt` file containing the standard 10X barcodes, which
    is available from 10X genomics and the LRSIM [GitHub repository](https://github.com/aquaskyline/LRSIM/).  
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock  --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/simulate-reads.smk')
    command.append('--configfile')
    command.append(f'{workflowdir}/config.yml')
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit(0)

    validate_input_by_ext(genome_hap1, "GENOME_HAP1", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])
    validate_input_by_ext(genome_hap2, "GENOME_HAP2", [".fasta", ".fa", ".fasta.gz", ".fa.gz"])

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fetch_file("simulate-reads.smk", f"{workflowdir}/")
    fetch_file("10xtoHaplotag.py", f"{workflowdir}/scripts/")
    fetch_file("LRSIMharpy.pl", f"{workflowdir}/scripts/")
    fetch_file("faFilter.pl", f"{workflowdir}/scripts/")
    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"genome_hap1: {genome_hap1}\n")
        config.write(f"genome_hap2: {genome_hap2}\n")
        config.write(f"output_directory: {output_dir}\n")
        if barcodes:
            config.write(f"barcodes: {barcodes}\n")
        config.write(f"outer_distance: {outer_distance}\n")
        config.write(f"distance_sd: {distance_sd}\n")
        config.write(f"read_pairs: {read_pairs}\n")
        config.write(f"molecule_length: {molecule_length}\n")
        config.write(f"partitions: {partitions}\n")
        config.write(f"molecules_per_partition: {molecules_per}\n")
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    onstart_text = f"Genome Haplotype 1: {os.path.basename(genome_hap1)}\n"
    onstart_text += f"Genome Haplotype 2: {os.path.basename(genome_hap2)}\n"
    onstart_text += f"Barcodes: {os.path.basename(barcodes)}\n" if barcodes else "Barcodes: 10X Default\n"
    onstart_text += f"Output Directory: {output_dir}/"
    print_onstart(onstart_text, "simulate reads")
    
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

