# MASHt - MASH toolkit
[//]: # (#TODO update readme mentioning the requirement to have term to group by as -1 in go_mart_to_go_csvs)

Python toolkit for the MASH tool and streamlined statistical analysis for Linux.

## Description

This toolkit serves to automate work that includes MASH [Ondov et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x). 

The blast module allows allows for multiprocessed automatized separation of transcripts/ reads by GO (gene ontology) terms. 

It also streamlines PCoA analysis on MASH distance matrices and performing ANOVA/MANOVA analysis on the results of the aforementioned separation to determine if there are statistically significant differences in gene expression levels between different groups of observations.

## Getting started

### Dependencies

- Python (tested on 3.10 and 3.11)
  - pandas
  - matplotlib
  - rpy2
  - statsmodels
  - skbio
- R
  - stats
  - broom
- mash (either in binary or installed in the environment)
- Linux OS (commands from the `stats` subcommand can also be run on Windows)

### Installation

No installation is necessary - the package should work right away.

Remember to install all dependencies.

## Executing the program
### Logical steps of analysis

A typical order of full analysis would be as follows:

1. BLAST step:
    1. downloading BioMart sequences and features (not necessary when using custom sequences and features)
    2. creating a BLAST database
    3. running BLAST
    4. splitting BLAST results by GOs
2. MASH step:
    1. creating MASH sketches
    2. calculating distances between sketches
3. stats step:
    1. performing PCoA analysis
    2. performing ANOVA/MANOVA analysis

The `--analyze_all` option for each respective module (step) allows for fully automated, multiprocessed analysis of specified files.

Since MASHt can take in arguments in a file preceded by `@` (e.g. `@args.txt`), it is recommended to use such files to make managing the arguments easier. Example files are provided in the `tests/` directory. 

<strong> NB. There should be no trailing spaces on any of the lines in the file! </strong> 


### mash module

- `--analyze_all` option allows for fully automated, multiprocessed analysis of masht.blaster `--split` results. MASH sketches files and triangle distance matrices will be created in all subdirectiories of selected file, e.g.:
    ```console
    foo@bar: python3 masht mash <path_to_dir> --analyze_all
    ```

- basic mash distance calculation (between the first and all of the rest of files in a specified directory):

    ```console
    foo@bar: python3 masht mash <path_to_dir> -d
    ```

    MASHt can also calculate distances between the first and the rest of specified files if `<path_to_dir>` is replaced with a `<path_to_file>`:

    ```console
    foo@bar: python3 masht mash <path_to_file> -d
    ```

    The `<path_to_file>` file should list filenames relative to current directory, e.g.:

    ``` txt
    tests/data/seqs/1.fastq
    tests/data/seqs/3.fastq
    tests/data/seqs/5.fastq
    ```

    NB - the `<path_to_file>` file has to have a .txt extension!
- basic mash sketching (of all files in the specified directory into one sketch):

    ```console
    foo@bar: python3 masht mash <path_to_dir> -s
    ```

- mash sketching, saving to specified folder and viewing information on these sketches:

    ```console
    foo@bar: python3 masht mash <path_to_dir> -s -i -o <output_dir>
    ```

- viewing info on all sketches in a folder:

    ```console
    foo@bar: python3 masht mash <path_to_dir> -i
    ```

- calculating distances between all sequences in a selected sketch (lower portion of the distance matrix) and printing it to console:

    ```console
    foo@bar: python3 masht mash <chosen_.msh_file> -t -v
    ```

- `--info` and `--triangle` options automatically work on the sketch file generated in the same command, e.g.:

    ```console
    foo@bar: python3 masht mash <dir_with_input_sequences> -s -t -i -o <output_dir>
    ```

    will create the sketch (.msh) file of all files within `dir_with_input_sequences/`, show information on them to the console and generate a `sketches_triangle.tsv` report file. The files will be stored in the `output_dir/`.
- parameters and options can be specified in a text file. Use `@<file_name>` to point to the file:

    ```console
    foo@bar: python3 masht mash @<file_name>
    ```

    `<file_name>` is relative to current directory, so:

    ```console
    foo@bar: python3 masht mash @test/args.txt
    ```

    will run masht mash with arguments provided in `args.txt` file within `./tests/` directory.
    `args.txt` should look something like this (note that there should be no trailing spaces on any of the lines!):

    ```text
    tests/data/seqs
    -s
    -o test_outputs
    -v
    ```

#### mash options (flags)

|option|long name|description|
|---|---|---|
|`-a`|`--analyze_all`|analyze .fasta files from all subdirectories of a given directory containing results of splitting blast results by GOs (blaster --split results)|
|`-b`|`--bounds`|calculate Mash error bounds of selected files|
|`-d`|`--distance`|calculate mash distances between the first and all of the rest of selected files|
|`-h`|`--help`|show help message|
|`-i`|`--info`|show information on selected sketch files|
|`-m`|`--mash`|run mash with specified params (point to a mash binary in the beginning)|
|`-o`|`--output_dir`|location of output directory (default: '.')|
|`-p`|`--paste`|paste multiple sketch files into a new one|
|`-s`|`--sketch`|create sketches of selected files|
|`-sc`|`--screen`|determine whether query sequences are within a sketch file|
|`-t`|`--triangle`|generate matrix of distances in a sketch|
|`-v`|`--verbose`|print more descriptions of performed actions to the console|

### stats module

- basic PCoA analysis:

    The simplest use case is to perform analysis on one triangle file generated by masht.mash `--triangle`:

    ```console
    foo@bar: python3 masht stats <path_to_triangle_file> -p -o <output_dir>
    ```

    One can also perform PCoA analysis on all files within `<path_to_dir>`:

    ```console
    foo@bar: python3 masht stats <path_to_dir> -p -o <output_dir>
    ```

    MASHt can also perform PCoA analysis on specified files if `<path_to_dir>` is replaced with a `<path_to_file>`:

    ```console
    foo@bar: python3 masht stats <path_to_file> -p
    ```

    The `<path_to_file>` file should list filenames relative to current directory, e.g.:

    ``` txt
    tests/data/triangle_one.tsv
    tests/data/triangle_two.tsv
    tests/data/triangle_three.tsv
    ```

    NB: the `<path_to_file>` file has to have a .txt extension!

  - number of dimensions is controlled with the `-n` option. Give an intiger to specify the number of dimensions or leave it empty to perform PCoA analysis on all dimensions:

    ```console
    foo@bar: python3 masht stats <path_to_triangle_file> -p -n 4
    ```

- ANOVA/MANOVA analysis:

    ```console
    foo@bar: python3 masht stats <path_to_pcoa_coords_file> -a -g <path_to_groups_file>
    ```

    **NB.** `<path_to_groups_file>` should be a tab-separated file!

    Names of observations in `<path_to_groups_file>` file should be the same as in `<path_to_pcoa_coords_file>` file. The file should look something like this:

    |Sample|Genotype|Temp|Time|
    |---|---|---|---|
    |BWOT1dR2|BW|OT|1d|
    |BWOT1dR1|BW|OT|1d|
    |BWOT1dR3|BW|OT|1d|
    |BWOT10dR1|BW|OT|10d|
    |BWOT10dR2|BW|OT|10d|
    |BWOT10dR3|BW|OT|10d|

  - `<path_to_pcoa_coords_file>` file is automatically inferred to be the results file of PCoA analysis is `-p` option was chosen, so there is no need to specify it explicitly. You still have to provide the `<path_to_groups_file>` file, though:

    ```console
    foo@bar: python3 masht stats <path_to_triangle_file> -p -a -g <path_to_groups_file>
    ```

    is equivalent to:

    ```console
    foo@bar: python3 masht stats <path_to_pcoa_coords_file> -a -g <path_to_groups_file>
    ```

    after the PCoA analysis is performed on `<path_to_triangle_file>` file.

    > NB: MANOVA analysis is performed using R via `rpy2` package to transpile code and objects (like dataframes) between R and Python. This makes the execution a bit slower, but testing showed that results of n-way MANOVA using `statsmodels` are not calculated correctly as of the time of writing this docs. Once the issue is resolved, the MANOVA analysis will be performed using statsmodels and execution time will improve.

- use the `--formula` option to specify the formula for ANOVA/MANOVA analysis. Only specify which exogenous (right hand side) variables to use and how to link them, e.g.: 'temp + genotype + time' - the endogenous (left hand side, PCs in this case) variables are inferred automatically based on other options, e.g. `--n`. <strong>By default a full model is created, e.g. one that includes all interactions between all the exogenous variables (with all exog variables with '*' between them)</strong>:

    ```console
    foo@bar: python3 masht stats <path_to_pcoa_coords_file> -a -g <path_to_groups_file> --formula 'Temp + Genotype + Time'
    ```

- parameters and options can be specified in a text file. Use `@<file_name>` to point to the file:

    ```console
    foo@bar: python3 masht stats @test/args.txt
    ```

    The file should look something like this (<strong>note that there should be no trailing spaces on any of the lines!</strong>):

    ```txt
    tests/anova_tests/aov_triangle.tsv
    -p
    -v
    -o
    test_outputs/
    -d
    1
    9
    -ma
    -g
    tests/anova_tests/test_groups_file.tsv
    ```

    this analysis will perform PCoA analysis on `tests/anova_tests/aov_triangle.tsv` file, draw a plot with PC1 and PC9 on it, perform MANOVA analysis on the same file and print more descriptions of performed actions to the console.

#### stats options (flags)

|option|long name|description|
|---|---|---|
||`--analyze_all`|perform PCoA and ANOVA/ MANOVA analysis on all files in the specified directory/ file|
|`-a`|`--anova`|perfom ANOVA on selected files. NB – a file with grouping has to be provided with the `-g` flag|
|`-amm`|`--anova_manova_mode`|select mode of ANOVA to perform. Should be either 'n' (to perform ANOVA on all parameters), an integer (for m-way ANOVA where first m columns from the groups_file will be selected) or 'repeat' for ANOVA with repeats. Defaults to 'n'
|`-d`|`--draw_plot`|draw PCoA plot when performing PCoA analysis. Two integers that signify which PCs to plot are required|
|`-f`|`--formula`|formula to use for ANOVA/MANOVA analysis. Specify only exogenous (right hand side, independent) variables and how to link them, e.g. 'temp + genotype + time'. The names must be consistent with columns in the `--groups_file` file.|
|`-g`|`--groups_file`|location of the .tsv file with groups for ANOVA/ MANOVA analysis. Required if `-a` or `-ma` was selected|
|`-m`|`--mode`|select which model to use for analysis when `--analyze_all` was chosen from [anova, manova]|
|`-ma`|`--manova`|perform MANOVA analysis on selected files. NB – a file with grouping has to be provided with the `-g` flag|
|`-n`|`-n_dimensions`|number of target dimensions for PCoA analysis|
|`-nt`|`--not_triangle`|signifies that the input file is NOT a triangle matrix (e.g., a distance matrix)|
|`-o`|`--output_dir`|location of output directory (default: '.')|
|`-p`|`--pcoa`|perform PCoA analysis on selected files|
|`-pc`|`--pcs`|number of PCs to analyse with ANOVA. Defaults to 4|
|`-ss`|`--ss_type`|type of sum of squares to use for ANOVA. Defaults to 2|
|`-v`|`--verbose`|print more descriptions of performed actions to the console|

### blaster module

All tasks related to BLAST (creation of index, blasting, splitting results by GO terms etc.) can be performed with one command. This way paths for input files for subsequent steps are inferred automatically. Example file with arguments would look like this:

```txt
--create_db
--db_fasta
../path/to/fasta_file_to_create_index_on (usually biomart sequences)
--name
name_of_the_index
--blast
--query
../path/to/file(s)_to_blast
--evalue
1e-49
--go_slim_list
../path/to/file(s)_with_go_terms
--split
../path/to/file(s)_with_sequences_to_split_by_go_terms (usually same as --query)
--output_dir
path/to/output_dir
--verbose
```

The same result can be achieved by using the ```--analyze_all``` option, template arguments file for which is provided in the `tests/` directory.


#### blaster options (flags)

|option|long name|description|
|---|---|---|
||`--analyze_all`|perform all BLAST-related tasks|
|`-b`|`--blast`|perform BLAST searches on files specified with selected `--query`|
|`-cdb`|`--create_db`|create BLAST database from the -db_fasta FASTA file|
|`-d`|`--download_biomart_files`|download GO Mart files from Ensembl Biomart based on data/*.xml queries|
|`-dbt`|`--db_type`|type of the BLAST database to create (nucl or prot). Default: nucl.|
|`-e`|`--evalue`|e-value threshold for BLAST search. Default: 10e-50|
|`-gsl`|`--go_slim_list`|create GO slim lists from provided GO Mart file|
|`-h`|`--help`|show this help message and exit
|`-n`|`--name`|name of the BLAST database to create or use|
|`-o`|`--output_dir`| output directory for the results|
|`-outfmt`|`--outfmt`|output format for BLAST results. Default: 6. Uses BLAST+ format codes.|
|`-q`|`--query`|query file to use for BLAST searches. Can either be a FASTA file or a file pointing to FASTA files (one per line) or a folder with FASTA files|
|`-s`|`--split`|split FASTA file provided here by GOs and BLAST results|
|`-v`|`--verbose`|verbose output|
||`--db_dir`|location of the BLAST database to use. Inferred automatically if `--create_db` is used|
||`--db_fasta`|location of the FASTA file to use for creating BLAST database|
||`--go`|location of the folder with go_list subfolders created by --go_slim_list. Inferred automatically if `--go_slim_list is used`|
||`--go_mart_feats`|path to file with GO features to use in --go_slim_list|
||`--in_blast_file`|location of the BLAST results file(s). Inferred automatically if `--blast` is used|
||`--n_jobs`|number of jobs to run in parallel for `--go_slim_list`. Default: 10|
||`--no_parse_seqids`|DO NOT parse SeqIDs in FASTA file when creating BLAST database|
||`--num_threads`|number of threads to use for BLAST search. Default: 4|


## Help

To access help, run the package with either -h or --help, e.g.:

```console
foo@bar: python3 masht -h
```

To print help for specific subcommand, specify it in your command:

``` console
foo@bar: python3 masht mash -h
```
