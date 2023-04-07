from masht.mash import _get_files, _print_error
import subprocess
import pathlib
import pandas as pd


def _read_fasta(fasta_file_path: str, prep: bool = False) -> dict[str, str]:
    """Read in fasta file and return a dict with the header as key and the sequence as value

    Args:
        fasta_file_path (str): path to fasta file

    Returns:
        dict[str, str]: dict with header as key and sequence as value
    """

    fasta_dict = {}

    if pathlib.Path(fasta_file_path).stat().st_size / 1e9 < 7:
        # faster way
        with open(fasta_file_path, 'r') as fasta_file:
            txt = fasta_file.read().split('>')[1:]
            for line in txt:
                line = line.split('\n')
                # no need for joining the sequence
                fasta_dict[line[0].split(' ')[0]] = line[1:]

    else:
        # memory efficient way
        with open(fasta_file_path, 'r') as fasta_file:
            curr_header = False
            for line in fasta_file:
                if line.startswith('>'):
                    if curr_header:
                        fasta_dict[header].append('')
                    header = line[1:].split(' ')[0]
                    fasta_dict[header] = []
                    curr_header = True
                else:
                    fasta_dict[header].append(line.rstrip())

    return fasta_dict


def split_blast_res_by_gos(blast_file_path: str or list[str], seqs_file_path: str, go_file_path: str or list[str], blast_outfmt: str = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', output_dir: str = '.', verbose: bool = 'False') -> None:
    """Split sequences from seqs_file based on go_file and blast_file. Puts resulting files in appropriate folders based on gene ontologies from the go_file.

        Args:
            blast_file_path (str or list[str]): path to blast file or list of paths to blast files. NB provide appropriate format for blast_outfmt as a string with qseqid, sseqed and pident at least!
            seqs_file_path (str): path to fasta file with sequences
            go_file_path (str or list[str]): path to file(s) with Gen Ontologies or list of paths to go files (from previous step)
            blast_outfmt (str, optional): blast output format. Defaults to 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' (standard outfmt 6 from BLAST).
            output_dir (str, optional): path to output directory. Defaults to '.'.
            verbose (bool, optional): whether to increase verbosity. Defaults to 'False'.
    """

    print('Splitting blast results by GOs...')

    # TODO remember to add the info that we require outfmt 6 to the docs

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

            blast_df.columns = blast_outfmt.split()

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


def blast_run(input_path: str, db: str, db_dir: str = '.', blast_type: str = 'blastn', evalue: float = 10e-50, num_threads: int = 4, outfmt: str = '6', output_dir: str = '.', verbose: bool = False) -> list[str]:
    # generate docs for this function
    """Run blast on a fasta file

    Args:
        input_path (str): path to fasta file
        db (str): name of blast database
        db_dir (str, optional): path to directory with blast database. Defaults to '.'.
        blast_type (str, optional): blast type. Defaults to 'blastn'.
        evalue (float, optional): evalue threshold. Defaults to 10e-50.
        num_threads (int, optional): number of threads. Defaults to 4.
        outfmt (str, optional): blast output format. Defaults to '6'. NB use a format with qseqid, sseqed and pident at least, since split_blast_res_by_gos requires those fields!
        output_dir (str, optional): path to output directory. Defaults to '.'.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.

    Returns:
        list[str]: list of paths to blast output files
    """

    if verbose:
        print('Running blast...')

    # assuming all inputs are in the input_dir folder
    in_files = _get_files(pathlib.Path(input_path))

    blast_files = []
    for file in in_files:
        proc = subprocess.run([blast_type, '-query', file, '-db', db, '-out', f'{file.stem}.blast', '-evalue', str(
            evalue), '-num_threads', str(num_threads), '-outfmt', outfmt], capture_output=True, cwd=pathlib.Path(db_dir))

        if proc.returncode != 0:
            _print_error(proc, 'blaster')
            return

        if verbose:
            print(proc.stdout.decode(), end='')

        blast_files.append(f'{db_dir}/{file.stem}.blast')

    return blast_files


def blast_create_index(input_file: str, name: str, db_type: str = 'nucl', parse_seqids: bool = True, verbose: bool = False) -> str:
    # generate docs for this function
    """Create blast index

    Args:
        input_file (str): path to input fasta file
        name (str): name of blast database
        db_type (str, optional): blast database type. Defaults to 'nucl'.
        parse_seqids (bool, optional): whether to parse seqids. Defaults to True.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.

    Returns:
        str: path to directory with blast database
    """

    print(f'Creating blast index for {input_file}...')

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

    return str(pathlib.Path(input_file).parent)  # blastdb location


def go_mart_to_go_slim_lists(go_file: str, output_dir: str) -> list[str]:
    """Split GO mart file to GO slim files in appropriate folders

    Args:
        go_file (str): path to GO mart file
        output_dir (str): path to output directory

    Returns:
        list[str]: list of paths to GO slim lists
    """

    print(f'Creating GO slim lists from {go_file}...')

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


# only for testing
if __name__ == '__main__':

    '''
    db_dir = blast_create_index(input_file='../testing/Hv_all_isoforms_sequence_mart_export.fasta',
                                name='test', db_type='nucl', parse_seqids=True, verbose=True)

    blast_files = blast_run(input_path='../testing/BWHT1dR3.fasta', db='test', db_dir='../testing')

    go_file = go_mart_to_go_slim_lists(
        go_file='../server_data/GO_skrypty/Hv_all_isoforms_GO_mart_export.txt', output_dir='../testing')
    '''

    split_blast_res_by_gos(blast_file_path='../testing/BWHT1dR3.blast', seqs_file_path='../testing/BWHT1dR3.fasta',
                           go_file_path='../testing/go_lists', output_dir='test_outputs')
