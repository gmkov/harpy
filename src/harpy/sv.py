from .harpymisc import getnames_err
import rich_click as click
import subprocess
import sys
import os

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-g', '--genome', type=click.Path(exists=True), required = True, metavar = "File Path", help = 'Genome assembly for variant calling')
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-p', '--populations', type=click.Path(exists = True), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def leviathan(genome, threads, directory, populations, extra_params, snakemake, quiet):
    """
    Call structural variants using LEVIATHAN
    
    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy extra --popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    samplenames = getnames_err(directory, '.bam')
    vcaller = "leviathan"
    if populations is not None:
        vcaller += "-pop"
    directory = directory.rstrip("/^")
    command = (f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/variants-{vcaller}.smk').split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    if populations is not None:
        # check that samplenames and populations line up
        with open(populations, "r") as f:
            p_list = [i.split()[0] for i in f.readlines()]
        missing_samples = [x for x in p_list if x not in samplenames]
        overlooked = [x for x in samplenames if x not in p_list]
        if len(missing_samples) > 0:
            print(f"\n\033[1;33mERROR:\033[00m There are {len(missing_samples)} samples included in \033[01m{populations}\033[00m that weren\'t found in \033[01m{directory}\033[00m. Terminating Harpy to avoid downstream errors. The samples causing this error are:", file = sys.stderr)
            print(", ".join(sorted(missing_samples)), file = sys.stderr)
            print(f"\n\033[1;34mSOLUTION:\033[00m Make sure the spelling of these samples is identical in \033[01m{directory}\033[00m and \033[01m{populations}\033[00m, or remove them from \033[01m{populations}\033[00m.\n", file = sys.stderr)
        if len(overlooked) > 0 and not quiet:
            print(f"\033[01mNOTICE:\033[00m There are {len(overlooked)} samples found in \033[01m{directory}\033[00m that weren\'t included in \033[01m{populations}\033[00m. This will not cause errors and can be ignored if it was deliberate. The samples are:", file = sys.stderr)
            print(", ".join(overlooked), file = sys.stderr)
        if len(missing_samples) > 0:        
            sys.exit(1)
        command.append(f"groupings={populations}")
    command.append(f"genomefile={genome}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True)
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with BAM alignments')
@click.option('-g', '--genome', required = True, type=click.Path(exists=True), metavar = "File Path", help = 'Genome assembly')
@click.option('-v', '--vcf', type=click.Path(exists=True), metavar = "File Path", help = 'Path to phased bcf/vcf file')
@click.option('-p', '--populations', type=click.Path(exists = True), metavar = "File Path", help = 'Tab-delimited file of sample<tab>population (optional)')
@click.option('-m', '--molecule-distance', default = 100000, show_default = True, type = int, metavar = "Integer", help = 'Base-pair distance delineating separate molecules')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional variant caller parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def naibr(genome, vcf, threads, directory, populations, molecule_distance, extra_params, snakemake, quiet):
    """
    Call structural variants using NAIBR
    
    NAIBR requires **phased** bam files as input. This appears as the `HP` or `PS` tags
    in alignment records. If your bam files do not have either of these phasing tags
    (e.g. BWA does not phase alignments), then provide a **phased** `--vcf` file such
     as that created by `harpy phase` and Harpy will use [whatshap haplotag](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
    to phase your input bam files prior to calling variants with NAIBR.

    Optionally specify `--populations` for population-pooled variant calling. 
    Use **harpy extra --popgroup** to create a sample grouping file to 
    use as input for `--populations`.
    """
    samplenames = getnames_err(directory, '.bam')
    vcaller = "naibr"
    if populations is not None:
        vcaller += "-pop"
    directory = directory.rstrip("/^")
    if vcf is not None:
        #TODO validate that the VCF is phased
        command = (f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/variants-{vcaller}-phase.smk').split()
    else:
        command = (f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/variants-{vcaller}.smk').split()
    
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"seq_directory={directory}")
    command.append(f"samplenames={samplenames}")
    if populations is not None:
        # check that samplenames and populations line up
        with open(populations, "r") as f:
            p_list = [i.split()[0] for i in f.readlines()]
        missing_samples = [x for x in p_list if x not in samplenames]
        overlooked = [x for x in samplenames if x not in p_list]
        if len(missing_samples) > 0:
            print(f"\n\033[1;33mERROR:\033[00m There are {len(missing_samples)} samples included in \033[01m{populations}\033[00m that weren\'t found in \033[01m{directory}\033[00m. Terminating Harpy to avoid downstream errors. The samples causing this error are:", file = sys.stderr)
            print(", ".join(sorted(missing_samples)), file = sys.stderr)
            print(f"\n\033[1;34mSOLUTION:\033[00m Make sure the spelling of these samples is identical in \033[01m{directory}\033[00m and \033[01m{populations}\033[00m, or remove them from \033[01m{populations}\033[00m.\n", file = sys.stderr)
        if len(overlooked) > 0 and not quiet:
            print(f"\033[01mNOTICE:\033[00m There are {len(overlooked)} samples found in \033[01m{directory}\033[00m that weren\'t included in \033[01m{populations}\033[00m. This will not cause errors and can be ignored if it was deliberate. The samples are:", file = sys.stderr)
            print(", ".join(overlooked), file = sys.stderr)
        if len(missing_samples) > 0:        
            sys.exit(1)
        command.append(f"groupings={populations}")
    command.append(f"molecule_distance={molecule_distance}")
    if vcf is not None:
        command.append(f"vcf={vcf}")
    if genome is not None:
        command.append(f"genomefile={genome}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)