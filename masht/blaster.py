from masht.mash import _get_files, _print_error
import subprocess
import pathlib
import pandas as pd

# TODO implement fasta preparation function that does the same as "awk '{print $1}' <fasta_file>"


def _read_fasta(fasta_file_path: str) -> dict[str, str]:
    """Read in fasta file and return a dict with the header as key and the sequence as value

    Args:
        fasta_file_path (str): path to fasta file

    Returns:
        dict[str, str]: dict with header as key and sequence as value
    """

    # faster way
    fasta_dict = {}
    with open(fasta_file_path, 'r') as fasta_file:
        txt = fasta_file.read().split('>')[1:]
        for line in txt:
            line = line.split('\n')
            # no need for joining the sequence
            fasta_dict[line[0]] = line[1:]
    return fasta_dict

    '''memory efficient way
    fasta_dict = {}
    with open(fasta_file_path, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                header = line[1:].rstrip()
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.rstrip()
    return fasta_dict
    '''


def split_blast_res_by_gos(blast_file_path: str or list[str], seqs_file_path: str, go_file_path: str or list[str], blast_outfmt: str = '6 qseqid sseqid pident evalue mismatch qstart sstart', output_dir: str = '.', verbose: bool = 'False') -> None:

    # TODO test this function
    '''
    blast_file_path: path to blast result file or list of paths to blast result files

    '''

    # TODO add option to use results from the previous steps as input for blast_file_path and go_file_path

    # get file lists to loop over
    blast_files = _get_files(pathlib.Path(blast_file_path)) if not isinstance(
        blast_file_path, list) else map(pathlib.Path, blast_file_path)
    go_files = _get_files(pathlib.Path(go_file_path)) if not isinstance(
        go_file_path, list) else map(pathlib.Path, go_file_path)
    seqs_files = _get_files(pathlib.Path(seqs_file_path))

    for go in go_files:
        # read in go file
        go_df = pd.read_csv(go)

        # create tun
        tun = (go_df['Gene stable ID'] + '|' +
               go_df['Transcript stable ID']).unique()

        # create output dir
        pathlib.Path(f'{output_dir}/{go.stem}').mkdir(
            parents=True, exist_ok=True)

        for blast_file in blast_files:
            blast_df = pd.read_csv(blast_file, sep='\t', header=None)
            # TODO make sure that outfmt in blast_run has correct columns (qseqid, pident) -> would be guaranteed if we just used the default 6 outfmt: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
            blast_df.columns = blast_outfmt.split()[1:]
            '''
            qseqid = blast_outfmt.split().index('qseqid')
            pident = blast_outfmt.split().index('pident')
            sseqid = blast_outfmt.split().index('sseqid')
            '''

            # select only the rows with highest e-value for each group based on qseqid
            mask = blast_df.groupby('qseqid')['pident'].transform(
                max) == blast_df['pident']

            # filter out rows in blast_df with highest pident for a group not in tun
            filtered_df = blast_df[mask & blast_df['sseqid'].isin(tun)]

            # get index of the corresponding sequence file to read it
            idx = [x.stem for x in seqs_files].index(blast_file.stem)
            seq_file = _read_fasta(seqs_files[idx])

            # write the corresponding sequences from seq_file to output file
            with open(f'{output_dir}/{go.stem}/filtered_{blast_file.stem}.fasta', 'w') as output_file:
                # get unique qseqids
                ids = filtered_df['qseqid'].unique()
                for id in ids:
                    output_file.write(">{0}\n{1}".format(
                        id, '\n'.join(seq_file[id])))


def blast_run(input_path: str, db: str, db_dir: str = '.', blast_type: str = 'blastn', evalue: float = 10e-50, num_threads: int = 1, outfmt: str = '6 qseqid sseqid pident evalue mismatch qstart sstart', output_dir: str = '.', verbose: bool = False) -> None:
    # assuming all inputs are in the input_dir folder
    in_files = _get_files(pathlib.Path(input_path))

    for file in in_files:
        proc = subprocess.run([blast_type, '-query', file, '-db', db, '-out', f'{file.stem}.blast', '-evalue', str(
            evalue), '-num_threads', str(num_threads), '-outfmt', outfmt], capture_output=True, cwd=pathlib.Path(db_dir))

        if proc.returncode != 0:
            _print_error(proc, 'blaster')
        return

    # TODO add return value of a list of paths to the blast files


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


split_blast_res_by_gos(blast_file_path='../testing/BWHT1dR3.blast', seqs_file_path='../testing/BWHT1dR3.fasta',
                       go_file_path='../testing/go_terms', output_dir='test_outputs')
