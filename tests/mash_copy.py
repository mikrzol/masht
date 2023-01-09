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


def run(args) -> None:
    files = _get_files(pathlib.Path(args['in_d']))

    for key, val in args.items():
        print(f'{key} -> {val}')

        # TODO do this:
        '''
        c. konkatenuj z resztą
        d. powtarzaj po wszystkich argumentach podanych przez usera
        e. odpal mash z tymi parametrami
        '''
