import subprocess
import pathlib


def _multiproc_task(subdirs: list[pathlib.Path]):
    sketch_path = sketch('bin/', data_path=subdirs,
                         output_path=subdirs, verbose=verb)
    triangle('bin/', data_path=pathlib.Path(sketch_path),
             output_path=subdirs, verbose=verb)


def analyze_all(go_dir: str, verbose: bool = False):
    """Analyze all files in subdirectories of go_dir (created with blaster.split_blast_res_by_gos). 

    Args:
        go_dir (str): path to directory created with blaster.split_blast_res_by_gos
    """
    import multiprocessing

    # assuming (filtered) .fasta files are to be analyzed
    subdirs = list(
        set(f.parent for f in pathlib.Path(go_dir).rglob('*.fasta')))

    # TODO add fastq support

    # workaround for passing verbose to _multiproc_task
    global verb
    verb = verbose

    with multiprocessing.Pool() as pool:
        pool.map(_multiproc_task, subdirs)


def _error_present(proc: subprocess.CompletedProcess, masht_subcommand: str) -> bool:
    """helper function for checking and printing errors

    Args:
        proc (subprocess.CompletedProcess): a completed process
        masht_subcommand (str): name of the masht subcommand

    Returns:
        bool: whether an error occurred
    """
    err = False
    if proc.returncode != 0:
        print(f'\nMASHt {masht_subcommand} encountered an error, see below:\n')
        print(proc.stderr.decode())
        err = True
    return err


def _get_files(data_path: pathlib.Path) -> list[pathlib.Path]:
    """helper function for getting files

    Args:
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files

    Returns:
        list[str]: list of files to process
    """
    files = []
    if data_path.is_file():
        if data_path.suffix == '.txt':
            with data_path.open('r') as in_f:
                files = [pathlib.Path(f) for f in in_f.read().strip().split()]
        else:
            files = [data_path]
    else:
        files = list(data_path.iterdir())

    return files


def _format_triangle_output(text: str) -> str:
    """Format the triangle command output for .tsv structure

    Args:
        text (str): mash triangle command output 

    Returns:
        str: formatted (.tsv style) output
    """
    import pandas as pd
    import numpy as np
    import re

    # read the text
    text = [line.split('\t') for line in text.strip().split('\n')][1:]
    # get sample names
    samples = [line[0] for line in text]
    # get values - lower triangle
    vals = np.zeros((len(samples), len(samples)))
    for i, line in enumerate(text):
        for j, el in enumerate(line[1:]):
            vals[i, j] = float(el)

    # get the full matrix
    full_mtx = np.triu(vals.T, 1) + vals
    # convert to pandas df
    full_df = pd.DataFrame(full_mtx, index=samples, columns=samples, dtype=str)
    # now get only the upper triangle
    triangle_df = full_df.where(np.triu(np.ones(full_df.shape)).astype(bool))
    triangle_df = triangle_df.stack().str.strip().reset_index()
    # name the columns
    triangle_df.columns = ['seq_A', 'seq_B', 'distance']
    # get rid of the diagonal
    triangle_df = triangle_df[triangle_df['seq_A'] != triangle_df['seq_B']]

    # convert to string
    triangle_df = triangle_df.to_string(index=False).strip()
    # clean up by converting multiple spaces to tabs
    triangle_df = re.sub(' +', '\t', triangle_df)
    # for some dumb reason sometimes lines start with space (changed to tabs above) even after .strip() - need to remove this manually
    return [el.lstrip('\t') for el in triangle_df.split('\n')]


def _loop_over_all_sketch_files(data_path: pathlib.Path, bin_path: pathlib.Path, mash_cmd: str, query: pathlib.Path = '', verbose: bool = False, output_path: str = '.') -> None:
    """helper funcion for looping over sketch files and applying chosen mash command

    Args:
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files
        bin_path (pathlib.Path): location of the mash binary files
        mash_cmd (str): mash command to perform on each file
        query (pathlib.Path): location of the query file (for screen command only)
        verbose (bool, optional): whether to print more info into terminal. Defaults to False.
        output_path (str): location of the folder to put the results into
    """
    # dict for storing names of commands that only print to the console
    console_only = {
        'info': 1,
        'screen': 1
    }

    # get files
    files = _get_files(data_path)

    for file in files:
        if file.suffix == '.msh':
            if verbose:
                print(f'============= {file}: =============')
            if query:
                query_files = _get_files(query)
                proc = subprocess.run([f'{bin_path}mash', mash_cmd,
                                       file, *query_files],
                                      capture_output=True)
                if _error_present(proc, mash_cmd):
                    return

            else:
                proc = subprocess.run([f'{bin_path}mash', mash_cmd,
                                       file],
                                      capture_output=True)
                if _error_present(proc, mash_cmd):
                    return

            if console_only.get(mash_cmd) or verbose:
                print(f"{proc.stdout.decode('utf-8')}", end='')

            if not console_only.get(mash_cmd):
                if mash_cmd == 'triangle':
                    with open(f'{output_path}/{file.name.split(".")[0]}_{mash_cmd}.tsv', 'w') as out_f:
                        out_f.write('\n'.join(_format_triangle_output(
                            proc.stdout.decode('utf-8'))))
                else:
                    with open(f'{output_path}/{file.name.split(".")[0]}_{mash_cmd}.txt', 'w') as out_f:
                        out_f.write(proc.stdout.decode('utf-8'))


def bounds(bin_path: str, data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> None:
    """calculate error bounds of selected sketch file[s]

    Args:
        bin_path (str): path that lead to a mash binary
        data_path (pathlib.Path): location of the sketch files
        output_path (str, optional): location of the report file. Defaults to '.'.
        verbose (bool, optional): whether to give more info in the console. Defaults to False.
    """
    print('\nError bounds of selected file[s]:')
    _loop_over_all_sketch_files(data_path=data_path, bin_path=bin_path,
                                mash_cmd='bounds', verbose=verbose, output_path=output_path)


def dist(bin_path: str, data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> None:
    """calculate mash distance between wanted files and save result to {output_path}/distances.tsv

    Args:
        bin_path (str): list of paths with mash binary
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files
        output_path (str, optional):  location of the directory to put the results.tsv file into. Defaults to '.'
        verbose (bool, optional): increase verbosity. Defaults to False.
    """
    print('\nCalculating mash distances...')
    files = _get_files(data_path)
    proc = subprocess.run([f'{bin_path}mash', 'dist',
                           *files],
                          capture_output=True
                          )
    if _error_present(proc, 'dist'):
        return

    if verbose:
        print(f"{proc.stdout.decode('utf-8')}", end='')

    with open(f'{output_path}/distances.tsv', 'w') as out_f:
        out_f.write('seq_A,seq_B,mash_dist,p_val,matching_hashes\n')
        out_f.write(proc.stdout.decode('utf-8'))
    print('\nMash distances calculated!')


def info(bin_path: str, data_path: str or pathlib.Path):
    """show info on selected sketch (.msh) files

    Args:
        bin_path (str): path that lead to a mash binary
        data_path (str or pathlib.Path): location of the sketch files
    """

    print('\nInformation on selected files:')

    _loop_over_all_sketch_files(
        data_path, bin_path, 'info')


def sketch(bin_path: str, data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> str:
    """
    apply mash sketch on given files if data_path leads to a single file or all files in a dir if it leads to a dir \
            and save result to {output_path}/sketches.msh

    Args:
        bin_path (str) : list of paths with mash binary
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files
        output_path (str, optional): location of the directory to put the sketches.msh file into. Defaults to '.'
        verbose: increase verbosity
    """
    if verbose:
        print('\nCreating mash sketches...')

    files = _get_files(data_path)

    proc = subprocess.run([f'{bin_path}mash', 'sketch',
                           *files,
                           '-o', f'{output_path}/sketches'], capture_output=True)

    if _error_present(proc, 'sketch'):
        return

    sketch_path = f'{output_path}/sketches.msh'

    if verbose:
        print(proc.stdout.decode())

    if verbose:
        print('\nMash sketches created!')

    return sketch_path


def paste(bin_path: str, data_path: pathlib.Path, file_name: str, output_path: str = '.') -> None:
    """paste one multiple sketch files into a new one

    Args:
        bin_path (str): path that lead to a mash binary
        data_path (pathlib.Path): location of the sketch files
        file_name (str): name of the new file
        output_path (str, optional): name of the new sketch file. Defaults to '.'.
    """
    if data_path.is_file():
        # txt file
        if data_path.suffix != '.msh':
            files = _get_files(data_path)

            proc = subprocess.run([f'{bin_path}mash', 'paste', f'{output_path}/{file_name}',
                                   *files,
                                   ],
                                  capture_output=True)
            if _error_present(proc, 'paste'):
                return
        # .msh file
        else:

            proc = subprocess.run([f'{bin_path}mash', 'paste', f'{output_path}/{file_name}',
                                   data_path,
                                   ],
                                  capture_output=True)
            if _error_present(proc, 'paste'):
                return
    # dir
    else:
        files = [file for file in data_path.iterdir() if file.suffix == '.msh']

        proc = subprocess.run([f'{bin_path}mash', 'paste', f'{output_path}/{file_name}',
                               *files,
                               ],
                              capture_output=True)
        if _error_present(proc, 'paste'):
            return


def screen(bin_path: str, data_path: pathlib.Path, query: pathlib.Path):
    """determine whether query files are within selected .msh files

    Args:
        bin_path (str): path to binary of mash
        data_path (pathlib.Path): input location
        query (pathlib.Path): location of query file[s]. Can be a single file, directory or file with relative paths to selected files.
    """
    print('\nScreening selected files:')

    _loop_over_all_sketch_files(
        data_path, bin_path, 'screen', query=query)


def triangle(bin_path: str, data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> None:
    """create matrix of distances between all of the sequences in a sketch to all others in this file

    Args:
        bin_path (str): path to binary of mash
        data_path (pathlib.Path): input location
        output_path (str, optional): location of the output directory. Defaults to '.'.
        verbose (bool, optional): whether to give more info on performed operations to the console. Defaults to False.
    """
    if verbose:
        print('\nMatrix of distances between sequences in selected files:')
    _loop_over_all_sketch_files(data_path=data_path, bin_path=bin_path,
                                mash_cmd='triangle', verbose=verbose, output_path=output_path)
