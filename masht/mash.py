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
        if data_path.suffix != '.msh':
            with data_path.open('r') as in_f:
                files = [pathlib.Path(f) for f in in_f.read().strip().split()]
        else:
            files = [data_path]
    else:
        files = list(data_path.iterdir())

    return files


def _loop_over_all_sketch_files(data_path: pathlib.Path, bin_paths: list[pathlib.Path], mash_cmd: str, verbose: bool = False, output_path: str = '.') -> None:
    """helper funcion for looping over sketch files and applying chosen mash command

    Args:
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files
        bin_paths (list[pathlib.Path]): location of the mash binary files (TESTING - remove this later)
        mash_func (str): mash command to perform on each file
        verbose (bool, optional): whether to print more info into terminal. Defaults to False.
    """
    # dict for storing names of commands that only print to the console
    console_only = {
        'info': 1
    }

    # get files
    files = _get_files(data_path)

    # TESTING - remove paths loop later
    for path in bin_paths:
        for file in files:
            if file.suffix == '.msh':
                print(f'============= {file}: =============')
                proc = subprocess.run([f'{path}mash', mash_cmd,
                                       file],
                                      capture_output=True)

                if console_only.get(mash_cmd) or verbose:
                    print(f"{proc.stdout.decode('utf-8')}", end='')

                if not console_only.get(mash_cmd):
                    with open(f'{output_path}/{file.name.split(".")[0]}_{mash_cmd}.txt', 'w') as out_f:
                        out_f.write(proc.stdout.decode('utf-8'))


def bounds(bin_paths: list[str], data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> None:
    print('Error bounds of selected file[s]:')
    _loop_over_all_sketch_files(data_path=data_path, bin_paths=bin_paths,
                                mash_cmd='bounds', verbose=verbose, output_path=output_path)


def dist(bin_paths: list[str], data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> None:
    """calculate mash distance on between wanted files and save result to {output_path}/distances.tsv

    Args:
        paths (list[str]) [testing var]: list of paths with mash binary
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files
        output_path (str, optional):  location of the directory to put the results.tsv file into. Defaults to '.'
        verbose (bool, optional): increase verbosity. Defaults to False.
    """
    print('\nCalculating mash distances...')
    files = _get_files(data_path)
    for path in bin_paths:
        proc = subprocess.run([f'{path}mash', 'dist',
                               *files],
                              capture_output=True
                              )
        if verbose:
            print(f"{proc.stdout.decode('utf-8')}", end='')

        with open(f'{output_path}/distances.tsv', 'w') as out_f:
            out_f.write('seq_A,seq_B,mash_dist,p_val,matching_hashes\n')
            out_f.write(proc.stdout.decode('utf-8'))
    print('Mash distances calculated!')


def sketch(bin_paths: list[str], data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> str:
    """
    apply mash sketch on given files if data_path leads to a single file or all files in a dir if it leads to a dir \
            and save result to {output_path}/sketches.msh

    Args:
        paths (list[str]) [testing var]: list of paths with mash binary
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files
        output_path (str, optional): location of the directory to put the sketches.msh file into. Defaults to '.'
        verbose: increase verbosity
    """
    print('\nCreating mash sketches...')
    # TESTING - remove types later
    types = ['old', 'new']
    files = _get_files(data_path)

    # TESTING - remove path loop later
    for path, type in zip(bin_paths, types):
        subprocess.run([f'{path}mash', 'sketch',
                        *files,
                        '-o', f'{type}_sketches'
                        ],
                       capture_output=not verbose)

    pathlib.Path(output_path).mkdir(
        parents=True, exist_ok=True)

    # move the *_sketches.msh files to a desired location (mash can only generate the sketch file to ./ ) --> TESTING (acually, it does do that)
    sketch_path = ''

    # TESTING - bin loop will be removed later
    for type, bin in zip(types, bin_paths):
        subprocess.run(
            ['mv', f'{type}_sketches.msh',
                f'{output_path}/{type}_sketches.msh'],
            capture_output=not verbose)
        sketch_path = f'{output_path}/{type}_sketches.msh'

    print('Mash sketches created!')

    return sketch_path


def info(bin_paths: list[str], data_path: str or pathlib.Path):
    """show info on selected sketch (.msh) files

    Args:
        bin_paths (list[str]): list of paths that lead to a mash binary (TESTING)
        data_path (str or pathlib.Path): location of the sketch files
    """

    # TESTING can remove this types of messages later
    print('Information on selected files:')
    if type(data_path) is str:  # i.e. sketch was called prior
        data_path = pathlib.Path(data_path)

    _loop_over_all_sketch_files(
        data_path, bin_paths, 'info')


def triangle(bin_paths: list[str], data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> None:
    print('Matrix of distances between sequences in selected files:')
    _loop_over_all_sketch_files(data_path=data_path, bin_paths=bin_paths,
                                mash_cmd='triangle', verbose=verbose, output_path=output_path)
