from mash import _get_files, _error_present
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
        while not (first_line := fasta_file.readline()).startswith('>'):
            first_line = fasta_file.readline()
        header = first_line[1:].split(' ')[0].rstrip()
        fasta_dict[header] = []

        for line in fasta_file:
            if line.startswith('>'):
                fasta_dict[header].append('')
                header = line[1:].split(' ')[0].rstrip()
                fasta_dict[header] = []
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

    if verbose:
        print('Splitting blast results by GOs...\n')

    # get file lists to loop over
    blast_files = _get_files(pathlib.Path(blast_file_path)) if not isinstance(
        blast_file_path, list) else list(map(pathlib.Path, blast_file_path))
    go_files = _get_files(pathlib.Path(go_file_path)) if not isinstance(
        go_file_path, list) else list(map(pathlib.Path, go_file_path))
    seqs_files = _get_files(pathlib.Path(seqs_file_path))

    for go in go_files:
        # read in go file
        go_df = pd.read_csv(go)

        # create tun
        tun = (go_df['Gene stable ID'] + '|' +
               go_df['Transcript stable ID']).unique()

        if verbose:
            print(f'Splitting {go.stem} file with {len(tun)} IDs...')

        # create output dir
        pathlib.Path(f'{output_dir}/{go.stem}').mkdir(
            parents=True, exist_ok=True)

        for blast_file in blast_files:
            if verbose:
                print(f'Processing {blast_file.stem} file...')

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

            # TODO this sometimes generates empty files - need to stop it then

            # get unique qseqids
            ids = filtered_df['qseqid'].unique()
            if not ids:
                continue

            # write the corresponding sequences from seq_file to output file
            with open(f'{output_dir}/{go.stem}/filtered_{blast_file.stem}.fasta', 'w') as output_file:
                for id in ids:
                    output_file.write(">{0}\n{1}".format(
                        id, '\n'.join(seq_file[id])))


def blast_run(input_path: str, db: str, db_dir: str = '.', blast_type: str = 'blastn', evalue: float = 10e-50, num_threads: int = 4, outfmt: str = '6', output_dir: str = '.', verbose: bool = False) -> list[str]:
    """Run blast on a fasta file

    Args:
        input_path (str): path to fasta file(s)
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

    # assuming all inputs are in the input_path folder
    in_files = _get_files(pathlib.Path(input_path))

    # TODO this can be multiprocessed

    blast_files = []
    for file in in_files:
        print(f'BLASTing {file}...')
        proc = subprocess.run([blast_type, '-query', file, '-db', db, '-out', f'{file.stem}.blast', '-evalue', str(
            evalue), '-num_threads', str(num_threads), '-outfmt', outfmt], capture_output=True, cwd=pathlib.Path(db_dir))

        if _error_present(proc, 'blast'):
            return

        if verbose:
            print(proc.stdout.decode(), end='')

        blast_files.append(f'{db_dir}/{file.stem}.blast')

    return blast_files


def blast_create_index(input_file: str, name: str, db_type: str = 'nucl', no_parse_seqids: bool = False, verbose: bool = False) -> str:
    # generate docs for this function
    """Create blast index

    Args:
        input_file (str): path to input fasta file
        name (str): name of blast database
        db_type (str, optional): blast database type. Defaults to 'nucl'.
        no_parse_seqids (bool, optional): whether to NOT parse seqids. Defaults to False.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.

    Returns:
        str: path to directory with blast database
    """

    print(f'Creating blast index for {input_file} ...')

    # create variables for subprocess.run
    parse_seqids = '-parse_seqids' if not no_parse_seqids else ''

    # run makeblastdb
    # TESTING - removed -parse_seqs temporarily
    proc = subprocess.run(['makeblastdb', '-in', pathlib.Path(input_file), '-dbtype',
                           db_type, '-title', name, '-out', f'{name}', '-blastdb_version', '5'], cwd=pathlib.Path(input_file).parent, capture_output=True)

    if verbose:
        print(proc.stdout.decode(), end='')

    if proc.returncode != 0:
        print('\nMASHt blaster encountered an error, see below:\n')
        print(proc.stderr.decode())
        return

    return str(pathlib.Path(input_file).parent)  # blastdb location


def go_mart_to_go_slim_lists(go_file: str, output_dir: str, n_jobs: int = 10) -> list[str]:
    """Split GO mart file to GO slim files in appropriate folders

    Args:
        go_file (str): path to GO mart file
        output_dir (str): path to output directory
        n_jobs (int): number of jobs to perform in parallel (defaults to 10)

    Returns:
        list[str]: list of paths to GO slim lists
    """
    from multiprocessing import Manager
    from joblib import Parallel, delayed

    def _mp_task(name: str, group: pd.DataFrame, output_dir: str, go_slim_list):
        group.to_csv(f'{output_dir}/go_csvs/{name}.csv', index=False)
        go_slim_list.append(f'{output_dir}/go_csvs/{name}.csv')

    print(f'Creating GO slim lists from {go_file}...')

    # create output dir
    pathlib.Path(f'{output_dir}/go_csvs').mkdir(
        parents=True, exist_ok=True)

    # read in go file
    go_df = pd.read_csv(go_file, sep='\t')

    # group by go slim terms (assuming last column is go slim term)
    grouped = go_df.groupby(go_df.columns[-1])

    with Manager() as manager:
        go_slim_list = manager.list()
        Parallel(n_jobs=n_jobs)(delayed(_mp_task)(name, group, output_dir, go_slim_list)
                                for name, group in grouped)

        return list(go_slim_list)


def query_biomart(output_dir: str, verbose: bool = False) -> dict:
    import urllib.request
    import pathlib

    if pathlib.Path('masht/data/biomart_feats_query.xml').is_file() and pathlib.Path('masht/data/biomart_seqs_query.xml').is_file():
        files = {'biomart_feats': 'masht/data/biomart_feats_query.xml',
                 'biomart_seqs': 'masht/data/biomart_seqs_query.xml'}

        for name, f_path in files.items():
            with open(f_path, 'r') as request:
                request = request.read().replace(
                    '\t', '').replace('\n', '').replace(' ', '%20')
                if verbose:
                    print(
                        f'fetching {name} from BioMart. This can take a while...')
                urllib.request.urlretrieve(
                    f'https://plants.ensembl.org/biomart/martservice?query={request}', f"{output_dir}/{name}.txt")
            '''
            # clean the files
            with open(f'{output_dir}/{name}.txt', 'r') as f_h:
                cleaned_text = f_h.read().replace(',', '&')
            '''

    else:
        print('Missing one or more data/*.xml files. Please ensure both biomart_feats_query.xml and data/biomart_seqs_query.xml are in the data/ directory.')

    return dict(zip(['feats', 'seqs'], [f'{output_dir}/{name}.txt' for name in ['biomart_feats', 'biomart_seqs']]))


# only for testing
# change from mash import to from masht.mash import to be able to run this part
if __name__ == '__main__':

    query_files = query_biomart(
        '/mnt/d/IGR_temp/new_blaster_test/', verbose=True)

    db_dir = blast_create_index(input_file=query_files['seqs'],
                                name='test', db_type='nucl', no_parse_seqids=True, verbose=True)

    blast_files = blast_run(
        input_path='/mnt/d/IGR_temp/blaster_test/inputs.txt', db='test', db_dir=db_dir, verbose=True)

    go_file = go_mart_to_go_slim_lists(
        go_file=query_files['feats'], output_dir='/mnt/d/IGR_temp/new_blaster_test')

    split_blast_res_by_gos(blast_file_path=blast_files, seqs_file_path='/mnt/d/IGR_temp/blaster_test/inputs.txt',
                           go_file_path=go_file, output_dir='/mnt/d/IGR_temp/new_blaster_test/go_lists_results/')
