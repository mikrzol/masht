# MASHt - MASH toolkit

Python toolkit for the MASH tool for Linux. 

## Description
This toolkit serves to automate work that includes MASH.


## Getting started
### Dependencies
- Python (tested on 3.10)
- mash (either in binary or installed in the environment)
- Linux OS


### Installation
No installation is necessary - the package should work right away. 


## Executing the program
### mash module
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
    ```
    tests/data/seqs
    -s
    -o test_outputs
    -v
    ```

## Options (flags):
- `-d`, `--distance`: calculate mash distances between the first and all of the rest of selected files
- `-s`, `--sketch`: create sketches of selected files
- `-o`, `--output_dir`: location of output directory (default: '.')
- `-v`, `--verbose`: add more descriptions of performed actions
- `-m`, `--mash`: run mash with specified params (point to a mash binary in the beginning)



## Help
To access help, run the package with either -h or --help, e.g.:
```console
foo@bar: python3 masht -h
```
To print help for specific subcommand, specify it in your command:
``` console
foo@bar: python3 masht mash -h
```
