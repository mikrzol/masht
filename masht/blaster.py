import pandas as pd
import pathlib
import subprocess
from masht.mash import _get_files, _print_error


def blast_run(input_path: str, output_dir: str, db_dir: str, db: str, blast_type: str = 'blastn', evalue: float = 10e-50, num_threads: int = 1, verbose: bool = False) -> None:
    # assuming all inputs are in the input_dir folder
    in_files = _get_files(pathlib.Path(input_path))

    for file in in_files:
        proc = subprocess.run([blast_type, '-query', file, '-db', db, '-out', f'{file.stem}.blast', '-evalue', str(
            evalue), '-num_threads', str(num_threads), '-outfmt', '6 qseqid sseqid pident evalue mismatch qstart sstart'], capture_output=True, cwd=pathlib.Path(db_dir))

        if proc.returncode != 0:
            _print_error(proc, 'blaster')
        return


def blast_create_index(input_file: str, name: str, db_type: str = 'nucl', parse_seqids: bool = True, verbose: bool = False) -> str:
    blastdb_location = str(pathlib.Path(input_file).parent)

    # create variables for subprocess.run
    parse_seqids = '-parse_seqids' if parse_seqids else ''

    # run makeblastdb
    proc = subprocess.run(['makeblastdb', '-in', pathlib.Path(input_file), '-dbtype',
                           db_type, '-title', name, '-out', f'{name}', parse_seqids], cwd=pathlib.Path(input_file).parent, capture_output=True)

    if verbose:
        print(proc.stdout.decode(), end='')

    if proc.returncode != 0:
        print('\nMASHt blaster encountered an error, see below:\n')
        print(proc.stderr.decode())
        return

    return blastdb_location


def go_mart_to_go_slim_lists(go_file: str, output_dir: str) -> list[str]:
    # create output dir
    pathlib.Path(f'{output_dir}/go_lists').mkdir(
        parents=True, exist_ok=True)

    # read in go file
    go_df = pd.read_csv(go_file)

    # group by go slim terms (assuming last column is go slim term)
    grouped = go_df.groupby(go_df.columns[-1])

    # save to appropriate csv files
    go_slim_list = []
    for name, group in grouped:
        group.to_csv(f'{output_dir}/go_lists/{name}.csv', index=False)
        go_slim_list.append(f'{output_dir}/go_lists/{name}.csv')

    return go_slim_list


def query_for_go_terms() -> None:
    # TODO implement querying for latest GO plant terms file from GO server
    pass
