from .helperfunctions import get_samples_from_fastq

import os
import re
import sys
import glob
import subprocess
import rich_click as click
from rich import print
from rich.panel import Panel

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-p', '--popgroup', required = False, type=click.Path(exists=True), metavar = "Input folder Path", help = 'Create generic sample-group file using existing sample file names (fq.gz or bam) in provided folder')
@click.option('-s', '--stitch-params', type=str, metavar = "Output file name", help = 'Create template STITCH parameter file')
@click.option('-h', '--hpc', type = click.Choice(["slurm", "sge"], case_sensitive = False), help = 'Create HPC scheduling profile')
def extra(popgroup, stitch_params, hpc):
    """
    Create various optional/necessary input files

    With this command you can generate a sample grouping file (for variant calling),
    a templace STITCH parameter file (for imputation), and a HPC profile for running
    Harpy on a cluster. You can use any combination of options at a time. 
    """
    if popgroup is not None:
        click.echo('\033[1m' + "<><> Sampling Grouping File <><>" + '\033[0m', file = sys.stderr, color = True)
        try:
            samplenames = getnames(popgroup, '.bam')
            click.echo("No bam files detected. Searching for fastq files.", file = sys.stderr)
        except:
            samplenames = get_samples_from_fastq(popgroup)

        click.echo(f"Samples detected in {popgroup}: " + str(len(samplenames)), file = sys.stderr)
        fout = "samples.groups"
        if os.path.exists("samples.groups"):
            overwrite = input("File \'samples.groups\' already exists, overwrite (no|yes)?  ").lower()
            if (overwrite == "no") or (overwrite == "n"):
                fout = input("Please suggest a different name for the output file: ")
            elif (overwrite == "yes") or (overwrite == "y"):
                fout = "samples.groups"
        with open(fout, "w") as file:
            for i in samplenames:
                file.write(i + '\tpop1\n') 
        click.echo('Created sample population grouping file: ' + fout + '\nPlease review it, as all samples have been grouped into a single population\n', file = sys.stderr, color = True)

    if stitch_params is not None:
        click.echo('\033[1m' + "<><> STITCH Parameter File <><>" + '\033[0m', file = sys.stderr)
        with open(stitch_params, "w") as file:
            file.write('model\tusebx\tbxlimit\tk\ts\tngen\ndiploid\tTRUE\t50000\t10\t5\t50\ndiploid\tTRUE\t50000\t10\t1\t50\ndiploid\tTRUE\t50000\t15\t10\t100')
        click.echo(f"Created example parameter file: {stitch_params}", file = sys.stderr)
        click.echo("Modify the model parameters as needed, but " + '\033[1m' + "DO NOT" + '\033[0m' + " add/remove columns", file = sys.stderr, color = True)

    if hpc is not None:
        click.echo('\033[1m' + "<><> HPC Profile <><>" + '\033[0m')
        subprocess.run(["hpc_profile.py", hpc])
