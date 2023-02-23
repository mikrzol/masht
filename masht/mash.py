import subprocess
import pathlib


def _get_files(data_path: pathlib.Path) -> list[str]:
    """helper function for getting files

    Args:
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files

    Returns:
        list[str]: list of files to process
    """
    files = []
    if data_path.is_file():
        # TESTING - temporary (?) solution for when file leads to .msh
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
    from itertools import combinations

    text = [line.split('\t') for line in text.strip().split('\n')]
    samples = [line[0] for line in text[1:]]
    values = [el for line in text[2:] for el in line[1:]]

    combs = [f'{line[0]}\t{line[1]}' for line in combinations(samples, r=2)]
    final = ['seq_A\tseq_B\tdistance'] + \
        ['\t'.join(el) for el in zip(combs, values)]

    return final


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

    # TESTING - remove paths loop later
    for file in files:
        if file.suffix == '.msh':
            print(f'============= {file}: =============')
            if query:
                query_files = _get_files(query)
                proc = subprocess.run([f'{bin_path}mash', mash_cmd,
                                       file, *query_files],
                                      capture_output=True)
                if proc.returncode != 0:
                    print(
                        f'MASH {mash_cmd} encountered an error, see below:\n')
                    print(proc.stderr.decode())
                    return

            else:
                proc = subprocess.run([f'{bin_path}mash', mash_cmd,
                                       file],
                                      capture_output=True)
                if proc.returncode != 0:
                    print(
                        f'MASH {mash_cmd} encountered an error, see below:\n')
                    print(proc.stderr.decode())
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
    """calculate mash distance on between wanted files and save result to {output_path}/distances.tsv

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
    if proc.returncode != 0:
        print('\nMASH dist encountered an error, see below:\n')
        print(proc.stderr.decode())
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

    # TESTING can remove this types of messages later
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
    print('\nCreating mash sketches...')

    files = _get_files(data_path)

    # TESTING - remove path loop later
    proc = subprocess.run([f'{bin_path}mash', 'sketch',
                           *files,
                           '-o', f'sketches'
                           ],
                          capture_output=True)

    if proc.returncode != 0:
        print(f'Encountered {proc.returncode} error, see below:\n')
        print(proc.stderr.decode())
        print('\nTerminating sketching...')
        return

    if verbose:
        print(proc.stdout.decode())

    pathlib.Path(output_path).mkdir(
        parents=True, exist_ok=True)

    # move the *_sketches.msh files to a desired location (mash can only generate the sketch file to ./ ) --> TESTING (acually, it does do that)
    sketch_path = ''

    subprocess.run(
        ['mv', f'sketches.msh',
            f'{output_path}/sketches.msh'],
        capture_output=not verbose)
    sketch_path = f'{output_path}/sketches.msh'

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

            # TESTING - remove paths loop later
            proc = subprocess.run([f'{bin_path}mash', 'paste', f'{output_path}/{file_name}',
                                   *files,
                                   ],
                                  capture_output=True)
            if proc.returncode != 0:
                print(f'MASH paste encountered an error, see below:\n')
                print(proc.stderr.decode())
                return
        # .msh file
        else:
            # TESTING - remove paths loop later

            proc = subprocess.run([f'{bin_path}mash', 'paste', f'{output_path}/{file_name}',
                                   data_path,
                                   ],
                                  capture_output=True)
            if proc.returncode != 0:
                print(f'MASH paste encountered an error, see below:\n')
                print(proc.stderr.decode())
                return
    # dir
    else:
        files = [file for file in data_path.iterdir() if file.suffix == '.msh']
        # TESTING - remove paths loop later

        proc = subprocess.run([f'{bin_path}mash', 'paste', f'{output_path}/{file_name}',
                               *files,
                               ],
                              capture_output=True)
        if proc.returncode != 0:
            print(f'MASH paste encountered an error, see below:\n')
            print(proc.stderr.decode())
            return


def screen(bin_path: str, data_path: pathlib.Path, query: pathlib.Path):
    """determine whether query files are within selected .msh files

    Args:
        bin_path (str): path to binary of mash
        data_path (pathlib.Path): input location
        query (pathlib.Path): location of query file[s]. Can be a single file, directory or file with relative paths to selected files.
    """
    # TESTING can remove this types of messages later
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
    print('\nMatrix of distances between sequences in selected files:')
    _loop_over_all_sketch_files(data_path=data_path, bin_path=bin_path,
                                mash_cmd='triangle', verbose=verbose, output_path=output_path)
