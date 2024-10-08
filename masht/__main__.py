import subprocess
import argparse
import pathlib

import mash
import stats
import blaster


def perform_blaster(args: argparse.ArgumentParser) -> None:
    # TODO use error logic to stop the program if an error occurs (instead of just ordering operations)
    """perform the blaster subcommand

    Args:
        args (argparse.ArgumentParser): args created in the main function
    """

    # either change this to '' or add default val to --db_dir arg = '' (or None)
    db_dir = '.'

    if args.analyze_all:
        blaster.analyze_all(args)
        return

    if args.download_biomart_files:
        try:
            biomart_files = blaster.query_biomart(
                output_dir=args.output_dir, verbose=args.verbose)
            args.db_fasta = biomart_files['seqs']
            args.go_mart_feats = biomart_files['feats']
        except Exception as e:
            print(f'Error encountered: {e}')
            return

    if args.create_db:
        try:
            db_dir = blaster.blast_create_index(input_file=args.db_fasta, name=args.name,
                                                db_type=args.db_type or 'nucl', no_parse_seqids=args.no_parse_seqids, output_dir=args.output_dir, verbose=args.verbose)
        except Exception as e:
            print(f'Error encountered: {e}')
            return

    blast_files = []
    if args.blast:
        try:
            blast_files = blaster.blast_run(
                input_path=args.query,
                db=args.name,
                db_dir=args.db_dir or db_dir,
                evalue=float(args.evalue) or 1e-49,
                num_threads=int(args.num_threads),
                outfmt=args.outfmt or '6',
                output_dir=args.output_dir,
                verbose=args.verbose)
        except Exception as e:
            print(f'Error encountered: {e}')
            return

    go_files = []
    if args.go_slim_list:
        try:
            go_files = blaster.go_mart_to_go_csvs(
                go_file=args.go_mart_feats, output_dir=args.output_dir)
        except Exception as e:
            print(f'Error encountered: {e}')
            return

    if args.split:
        try:
            blaster.split_blast_to_fastas(
                blast_file_path=blast_files or args.in_blast_file,
                go_file_path=go_files or args.go,
                seqs_file_path=args.split,
                blast_outfmt=args.outfmt or 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
                output_dir=args.output_dir,
                verbose=args.verbose)
        except Exception as e:
            print(f'Error encountered: {e}')
            return


def perform_stats(args: argparse.ArgumentParser) -> None:
    # TODO use error logic to stop the program if an error occurs (instead of just ordering operations)
    """perform the stats subcommand

    Args:
        args (argparse.ArgumentParser): args created in the main function
        bin_path (str): directory with the mash binary (obsolete)
        data_path (pathlib.Path): location of the dir/ path with files to work on
    """

    data_path = pathlib.Path(args.in_d)

    # TESTING - remove paths later
    bin_path = 'bin/'

    if args.n_dimensions:
        args.n_dimensions = int(args.n_dimensions)
    if args.pc_number:
        args.pc_number = int(args.pc_number)

    args.ss_type = int(args.ss_type)

    if args.analyze_all:
        stats.analyze_all(data_path=data_path,
                          mode=args.mode,
                          groups_file=args.groups_file,
                          output_dir=pathlib.Path(args.output_dir),
                          anova_manova_mode=args.anova_manova_mode,
                          formula=args.formula,
                          pcs=args.pc_number,
                          verbose=args.verbose,
                          plot=args.draw_plot,
                          ss_type=args.ss_type,
                          triangle=args.not_triangle,
                          n_dim=args.n_dimensions)
        return

    pcoa_path = ''
    if args.pcoa:
        pcoa_path = stats.pcoa(data_path=data_path, output_dir=args.output_dir, n_dim=args.n_dimensions,
                               plot=args.draw_plot, triangle=args.not_triangle, verbose=args.verbose)
    if pcoa_path:
        pcoa_path = pathlib.Path(pcoa_path)

    if args.anova:
        stats.anova(data_path=pcoa_path or data_path, groups_file=args.groups_file,
                    anova_manova_mode=args.anova_manova_mode, formula=args.formula, output_dir=args.output_dir, pcs=args.pc_number, ss_type=args.ss_type, triangle=args.not_triangle, verbose=args.verbose)

    if args.manova:
        stats.manova(data_path=pcoa_path or data_path,
                     groups_file=args.groups_file, anova_manova_mode=args.anova_manova_mode, output_dir=args.output_dir, pcs=args.pc_number, verbose=args.verbose)


def perform_mash(args: argparse.ArgumentParser) -> None:
    # TODO use error logic to stop the program if an error occurs (instead of just ordering operations)
    """perform mash command

    Args:
        args (argparse.ArgumentParser): args created in the main function
        bin_path (str): directory with the mash binary
        data_path (pathlib.Path): location of the dir/ path with files to work on
    """

    data_path = pathlib.Path(args.in_d)

    # TESTING - remove paths later
    bin_path = 'bin/'

    # ORDER MATTERS
    # sketch
    sketch_path = ''
    if args.sketch:
        sketch_path = mash.sketch(bin_path=bin_path,
                                  data_path=data_path, output_path=args.output_dir,
                                  verbose=args.verbose)
    if sketch_path:
        sketch_path = pathlib.Path(sketch_path)
    # info
    if args.info:
        mash.info(bin_path=bin_path, data_path=sketch_path or data_path)

    # bounds
    if args.bounds:
        mash.bounds(bin_path=bin_path, data_path=sketch_path or data_path,
                    output_path=args.output_dir, verbose=args.verbose)
    # dist
    if args.distance:
        mash.dist(bin_path=bin_path, data_path=data_path,
                  output_path=args.output_dir, verbose=args.verbose)

    # triangle
    if args.triangle:
        mash.triangle(bin_path=bin_path, data_path=sketch_path or data_path,
                      output_path=args.output_dir, verbose=args.verbose)

    # paste
    if args.paste:
        mash.paste(bin_path=bin_path, data_path=data_path,
                   output_path=args.output_dir, file_name=args.paste)

    # screen
    if args.screen:
        mash.screen(bin_path=bin_path,
                    data_path=data_path, query=pathlib.Path(args.screen))

    # mash (fully custom)
    if args.mash:
        subprocess.run(args.mash.split())

    if args.analyze_all:
        mash.analyze_all(go_dir=data_path, verbose=args.verbose)


def main():
    global_parser = argparse.ArgumentParser(prog='MASHt',
                                            description='MASH distance toolkit',
                                            fromfile_prefix_chars='@')

    subparsers = global_parser.add_subparsers(
        title='subcommands', help='available subcommands', dest='subparser_name')

    # BLASTER SUBCOMMAND
    blaster_parser = subparsers.add_parser(
        'blaster', help='use the blaster module')

    # BLASTER - SHARED ARGUMENTS
    blaster_parser.add_argument(
        '--analyze_all', action='store_true', help='run all the analysis steps in one go. Most other options can still be used to customize the analysis.')

    blaster_parser.add_argument(
        '-o', '--output_dir', help='output directory for the results')

    blaster_parser.add_argument(
        '-outfmt', '--outfmt', help='output format for BLAST results. Default: 6. Uses BLAST+ format codes.')

    blaster_parser.add_argument(
        '-v', '--verbose', action='store_true', help='verbose output')

    blaster_parser.add_argument(
        '-n', '--name', help='name of the BLAST database to create or use')

    # BLASTER - QUERYING BIOMART FOR FEATS AND SEQS
    blaster_parser.add_argument('-d', '--download_biomart_files', action='store_true',
                                help='whether to fetch feature and sequences specified in data/*.xml files from BioMart.')

    # BLASTER - CREATING BLAST DATABASE
    blaster_parser.add_argument(
        '-cdb', '--create_db', action='store_true', help='create BLAST database from the --db_fasta FASTA file')

    blaster_parser.add_argument(
        '-dbt', '--db_type', help='type of the BLAST database to create (nucl or prot). Default: nucl')

    blaster_parser.add_argument('--no_parse_seqids', action='store_true',
                                help='DO NOT parse SeqIDs in FASTA file when creating BLAST database')

    blaster_parser.add_argument('--db_fasta',
                                help='location of the FASTA file to use for creating BLAST database')

    # BLASTING
    blaster_parser.add_argument(
        '-b', '--blast', action='store_true', help='perform BLAST searches on selected files')

    blaster_parser.add_argument(
        '--db_dir', help='location of the BLAST database to use. Inferred automatically if --create_db is used')

    blaster_parser.add_argument(
        '-q', '--query', help='query file to use for BLAST searches. Can either be a FASTA file, a file pointing to FASTA files (one per line) or a folder with FASTA files')

    blaster_parser.add_argument('-e', '--evalue', type=float, default=1e-49,
                                help='e-value threshold for BLAST search. Default: 1e-49')

    blaster_parser.add_argument('--num_threads', type=int, default=4,
                                help='number of threads to use for BLAST search. Default: 4')

    # CREATING GO SLIM LISTS FROM GO MART FILE
    blaster_parser.add_argument(
        '-gsl', '--go_slim_list', action='store_true', help='create GO slim lists from provided GO Mart file')

    blaster_parser.add_argument(
        '--go_mart_feats', help='path to file with GO features to use in --go_slim_list')

    blaster_parser.add_argument(
        '--n_jobs', default=10, help='number of jobs to perform in parallel for --go_slim_list')

    # SPLITTING FASTA FILE BY GOs AND BLAST RESULTS
    blaster_parser.add_argument(
        '-s', '--split', help='split FASTA file provided here by GOs and BLAST results')

    blaster_parser.add_argument(
        '--go', help='location of the folder with go_list subfolders created by --go_slim_list. Inferred automatically if --go_slim_list is used')

    blaster_parser.add_argument(
        '--in_blast_file', help='location of the BLAST results file(s). Inferred automatically if --blast is used')

    # set function to perform when calling the command
    blaster_parser.set_defaults(func=perform_blaster)

    # STATS SUBCOMMAND
    stats_parser = subparsers.add_parser('stats', help='use the stats module')

    '''
    stats_parser.add_argument(
        '--input', type=argparse.FileType('r'), help='input file')
    '''

    stats_parser.add_argument('in_d',
                              help='location of 1) the folder with FASTQ or FASTA \
                           files or 2) the file with names of selected files \
                               (names are relative to package location)')
    stats_parser.add_argument(
        '-a', '--anova', action='store_true', help='perfom ANOVA on selected files')
    stats_parser.add_argument('--analyze_all', action='store_true',
                              help='perform full analysis on all files in the data directory and subdirectories')
    stats_parser.add_argument(
        '-amm', '--anova_manova_mode', default='n', help='select mode of ANOVA to perform. Should be either \'n\' (to perform ANOVA on all parameters), an integer (for m-way ANOVA where first m columns from the groups_file will be selected) or \'repeat\' for ANOVA with repeats. Defaults to \'n\'')
    stats_parser.add_argument(
        '-d', '--draw_plot', nargs=2, help='draw PCoA plot for chosen PCs. Two intigers required.')
    stats_parser.add_argument('-f', '--formula', default=None,
                              help='formula to use for ANOVA or MANOVA. Required if -a or -ma was selected')
    stats_parser.add_argument(
        '-g', '--groups_file', help='location of the file containing information on grouping for ANOVA or MANOVA. Required if -a or -ma was selected')
    stats_parser.add_argument('-ma', '--manova', action='store_true',
                              help='perfom MANOVA analysis on selected files')
    stats_parser.add_argument('-m', '--mode', choices=[
                              'anova', 'manova'], help='select which model to use for analysis when `--analyze_all` was chosen.')
    stats_parser.add_argument('-n', '--n_dimensions', default=None,
                              help='number of target dimensions for PCoA analysis')
    stats_parser.add_argument('-nt', '--not_triangle', action='store_false',
                              help='signifies that the input file is a NOT a triangle matrix')
    stats_parser.add_argument('-o', '--output_dir', default='./',
                              help='location of the output directory (default: ".")')
    stats_parser.add_argument('-p', '--pcoa', action='store_true',
                              help='perform PCoA analysis and create results files')
    stats_parser.add_argument('-pc', '--pc_number', default=4,
                              help='Number of PCs to analyse with ANOVA. Defaults to 4.')
    stats_parser.add_argument(
        '-ss', '--ss_type', choices=['1', '2', '3'], default='2', help='Type of sum of squares for ANOVA.')
    stats_parser.add_argument('-v', '--verbose', action='store_true',
                              help='add more descriptions of performed actions')

    # set function to perform when calling the command
    stats_parser.set_defaults(func=perform_stats)

    # MASH SUBCOMMAND
    mash_parser = subparsers.add_parser('mash', help='use the mash module')
    basic = mash_parser.add_argument_group('basic use')
    basic.add_argument('in_d',
                       help='location of 1) the folder with FASTQ or FASTA \
                           files or 2) the file with names of selected files \
                               (names are relative to package location)')
    detailed = mash_parser.add_argument_group('detailed use')
    detailed.add_argument(
        '-a', '--analyze_all', action='store_true', help='analyze files from all subdirectories of a given directory containing results of splitting blast results by GOs (blaster --split results)')
    detailed.add_argument('-b', '--bounds', action='store_true',
                          help='show Mash error bounds of selected files')
    detailed.add_argument('-d', '--distance', action='store_true',
                          help='calculate distance between selected files')
    detailed.add_argument('-i', '--info', action='store_true',
                          help='show information on selected sketch files')
    detailed.add_argument('-m', '--mash', help='run mash with specified params',
                          type=str)
    detailed.add_argument('-o', '--output_dir', default='./',
                          help='location of the output directory (default: ".")')
    detailed.add_argument('-p', '--paste',
                          help='paste multiple sketch files into a new one')
    detailed.add_argument('-s', '--sketch', action='store_true',
                          help='generate mash sketches of selected files')
    detailed.add_argument('-sc', '--screen',
                          help='determine whether query sequences are within a sketch file')
    detailed.add_argument('-t', '--triangle', action='store_true',
                          help='generate matrix of distances in a sketch')
    detailed.add_argument('-v', '--verbose', action='store_true',
                          help='add more descriptions of performed actions')

    # set function to perform when calling the command
    mash_parser.set_defaults(func=perform_mash)

    # parse args
    '''
    # Parse command line arguments once to extract path to input file
    temp_args, remaining_args = global_parser.parse_known_args()
    if temp_args.input:
        with temp_args.input as f:
            contents = f.read()
        # Split contents on whitespace delimiters
        inputs = contents.split()
        # Append inputs to remaining command line arguments
        remaining_args.extend(inputs)
    # Parse command line arguments again with input file contents added
    args = global_parser.parse_args(remaining_args)
    '''

    args = global_parser.parse_args()
    if args.subparser_name == 'stats':
        if (bool(args.manova) or bool(args.anova) or bool(args.analyze_all)) ^ bool(args.groups_file):
            global_parser.error(
                '--anova or --manova or --analyze_all and --groups_file must be given together!')
            return

    if args.subparser_name == 'blaster':
        if (not (bool(args.download_biomart_files) or bool(args.analyze_all)) and (bool(args.go_slim_list) ^ bool(args.go_mart_feats))):
            global_parser.error(
                'either --download_biomart_files or --analyze_all must be used or --go_slim_list and --go_mart_feats must be given together!')
            return

    # VARIABLES
    # check if output dir exists and create one if necessary
    args.output_dir = args.output_dir.strip()
    pathlib.Path(f'{args.output_dir}').mkdir(
        parents=True, exist_ok=True)

    # MASH
    # perform_mash(args=args, bin_path=bin_path, data_path=data_path)

    # STATS
    # perform_stats(args=args, bin_path=bin_path, data_path=data_path)

    # perform appropriate function
    args.func(args)


if __name__ == '__main__':
    main()
