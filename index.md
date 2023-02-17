---
label: Home
icon: home
---
![](static/logo.png)

Harpy is a haplotagging data processing pipeline for Linux-based systems. It uses all the 
magic of [Snakemake](https://snakemake.readthedocs.io/en/stable/) under the hood to handle 
the worklfow decision-making, but as a user, you just interact with it like a normal command-line 
program. Harpy uses both well known and niche programs to take raw haplotagging sequences and process
them to become called SNP genotypes (or haplotypes). Most of the settings are pre-configured and the settings you
can modify are done at the command line. There aren't too many, which should make things a little simpler. 

## Modules
Harpy is modular, meaning you can use different parts of it independent from each other. Need to only align reads?
Great! Only want to call variants? Awesome! All modules are called by `harpy <module>`. For example, use `harpy align` to align reads.

| Module     | Description                                   |
|:-----------|:----------------------------------------------|
| `extra`    | Create various associated or necessary files  |
| `trim`     | Remove adapters and quality trim sequences    |
| `align`    | Align sample sequences to a reference genome  |
| `variants` | Call variants from sample alignments          |
| `impute`   | Impute genotypes using variants and sequences |
| `phase`    | Phase SNPs into haplotypes                    |


## Using Harpy
You can call `harpy` without any arguments (or with `--help`) to print the docstring to your terminal. You can likewise call any of the modules with `--help` (e.g. `harpy align --help`) to see their usage.
``` harpy --help                                                      
 Usage: harpy [OPTIONS] COMMAND [ARGS]...                     
                                                              
               Haplotagging Research Pipeline (HARPY)               
                            version: 0.1                            
                                                                    
 The pipeline trims reads, maps sequences, calls variants, imputes  
 genotypes, and phases haplotypes of Haplotagging data.             
                                                                    
 trim 🡒 align 🡒 variants 🡒 impute 🡒 phase                           
                                                                    
 Documentation: https://pdimens.github.io/HARPY/                    
                                                                    
╭─ Options ────────────────────────────────────────────────────────╮
│ --help      Show this message and exit.                          │
╰──────────────────────────────────────────────────────────────────╯
╭─ Commands ───────────────────────────────────────────────────────╮
│ align     Align sample sequences to a reference genome           │
│ extra     Create various associated/necessary files              │
│ impute    Impute genotypes using variants and sequences          │
│ phase     Phase SNPs into haplotypes                             │
│ trim      Remove adapters and quality trim sequences             │
│ variants  Call variants from sample alignments                   │
╰──────────────────────────────────────────────────────────────────╯
```