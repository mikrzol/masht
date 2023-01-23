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
    ```
    tests/data/seqs
    -s
    -o test_outputs
    -v
    ```

#### mash options (flags):
|option|long name|description|
|---|---|---|
|`-b`|`--bounds`|calculate Mash error bounds of selected files|
|`-d`|`--distance`|calculate mash distances between the first and all of the rest of selected files|
|`-i`|`--info`|show information on selected sketch files|
|`-m`|`--mash`|run mash with specified params (point to a mash binary in the beginning)|
|`-o`|`--output_dir`|location of output directory (default: '.')|
|`-s`|`--sketch`|create sketches of selected files|
|`-sc`|`--screen`|determine whether query sequences are within a sketch file|
|`-t`|`--triangle`|generate matrix of distances in a sketch|
|`-v`|`--verbose`|print more descriptions of performed actions to the console|


## Help
To access help, run the package with either -h or --help, e.g.:
```console
foo@bar: python3 masht -h
```
To print help for specific subcommand, specify it in your command:
``` console
foo@bar: python3 masht mash -h
```
