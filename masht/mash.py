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
        with data_path.open('r') as in_f:
            files = in_f.read().strip().split()
    else:
        files = list(data_path.iterdir())

    return files


def dist(bin_paths: list[str], data_path: pathlib.Path, output_path: str = '.', verbose: bool = False) -> None:
    """calculate mash distance on between wanted files and save result to {output_path}/distances.csv

    Args:
        paths (list[str]) [testing var]: list of paths with mash binary
        data_path (pathlib.Path): location of the file with selected files listed or the dir with wanted files
        output_path (str, optional):  location of the directory to put the results.csv file into. Defaults to '.'
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

        with open(f'{output_path}/distances.csv', 'w') as out_f:
            out_f.write('seq_A,seq_B,mash_dist,p_val,matching_hashes\n')
            out_f.write(proc.stdout.decode('utf-8').replace('\t', ','))
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

    # move the *_sketches.msh files to a desired location (mash can only generate the sketch file to ./ )
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
    if type(data_path) is str:  # i.e. sketch was called prior
        data_path = pathlib.Path(data_path)

    # TESTING - remove path loop later
    types = ['old', 'new']
    for path, typ in zip(bin_paths, types):
        subprocess.run([f'{path}mash', 'info',
                        data_path],
                       capture_output=False)
