---
label: Common Options
icon: list-unordered
order: 4
---

# :icon-list-unordered: Common Harpy Options
## Input Arguments
Each of the main Harpy modules (e.g. `qc` or `phase`) follows the format of
```bash
harpy module options arguments
```
where `module` is something like `impute` or `snp mpileup` and `options` are the runtime parameters,
which can include things like an input `--vcf` file, `--molecule-distance`, etc. After the options
is where you provide the input files/directories without flags and following standard BASH expansion
rules (e.g. wildcards). You can mix and match entire directories, individual files, and wildcard expansions.
In most cases, you can provide an unlimited amount of input arguments, which Harpy will parse and symlink
into the `*/workflow/input` folder, leaving the original files unmodified. In practice, that can look like:
```bash
harpy align bwa -t 5 -g genome.fasta data/pop1 data/pop2/trimmed*gz data/pop3/sample{1,2}* data/pop4/sample{2..5}*gz 
```
!!!info not recursive
Keep in mind that Harpy will not recursively scan input directories for files. If you provide `data/` as an input,
Harpy will search for fastq/bam files in `data/` and not in any subdirectories within `data/`. This is done deliberately
to avoid unexpected behavior.
!!!

!!!warning clashing names
Harpy will symlink just the file names into `workflow/input` regardless of their origin,
meaning that files in different directories that have the same name (ignoring extensions) will
clash. As an example, both `folderA/sample001.bam` and `folderB/sample001.bam` will become symlinked
as `workflow/input/sample001.bam`, with one symlink overwriting the other, leaving you with one missing
sample. During parsing, Harpy will inform you of naming clashes and terminate to protect you against
this behavior. 
!!!

## Common command-line options
Every Harpy module has a series of configuration parameters. These are arguments you need to input
to configure the module to run on your data, such as the directory with the reads/alignments,
the genome assembly, etc. All main modules (e.g. `qc`) also share a series of common runtime
parameters that don't impact the results of the module, but instead control the speed/verbosity/etc.
of calling the module. These runtime parameters are listed in the modules' help strings and can be 
configured using these arguments:

| argument        | short name | type    | default | required | description                                                                       |
|:--------------- |:----------:|:------- |:-------:|:--------:|:--------------------------------------------------------------------------------- |
| `--output-dir`  | `-o`       | string  | varies  | no       | Name of output directory                                                          |
| `--threads`     | `-t`       | integer | 4       | no       | Number of threads to use                                                          |
| `--skipreports` | `-r`       | toggle  |         | no       | Skip the processing and generation of HTML reports in a workflow                  |
| `--snakemake`   | `-s`       | string  |         | no       | Additional [Snakemake](snakemake/#adding-snakamake-parameters) options, in quotes |
| `--quiet`       | `-q`       | toggle  |         | no       | Supressing Snakemake printing to console                                          |
| `--help`        |            |         |         |          | Show the module docstring                                                         |

As as example, you could call the `harpy align` module and specify 20 threads with no output to console:

```bash
harpy align bwa --threads 20 --quiet samples/trimmedreads

# identical to #

harpy align bwa -t 20 -q samples/trimmedreads
```
---

## The `workflow` folder
When you run one of the main Harpy modules, the output directory will contain a `workflow` folder. This folder is
both necessary for the module to run and is very useful to understand what the module did, be it for your own
understanding or as a point of reference when writing the Methods within a manuscript. The presence of the folder
and the contents therein also allow you to rerun the workflow manually. The `workflow` folder may contain the following:

| item | contents | utility |
|:-----|:---------|:--------|
|`*.smk`               | Snakefile with the full recipe of the workflow | useful for understanding the workflow |
| `config.yml`         | Configuration file generated from command-line arguments and consumed by the Snakefile | useful for bookkeeping | 
| `input/`             | Symlinks to all of the provided input files with standardized extensions |
| `report/*.Rmd`       | RMarkdown files used to generate the fancy reports | useful to understand math behind plots/tables or borrow code from |
| `*.workflow.summary` | Plain-text overview of the important parts of the workflow | useful for bookkeeping and writing Methods |

---

## The `Genome` folder
You will notice that many of the workflows will create a `Genome` folder in the working 
directory. This folder is to make it easier for Harpy to store the genome and the associated
indexing/etc. files. Your input genome will be symlinked into that directory (not copied), but
all the other files (`.fai`, `.bwt`, `.bed`, etc.) will be created in that directory.