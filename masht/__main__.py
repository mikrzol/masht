import subprocess
import argparse
import pathlib

import mash
import stats
import blaster


def perform_blaster(args: argparse.ArgumentParser) -> None:
    """perform the blaster subcommand

    Args:
        args (argparse.ArgumentParser): args created in the main function
    """

    db_dir = '.'
    if args.create_db:
        db_dir = blaster.blast_create_index(input_file=args.db_fasta, name=args.name,
                                            db_type=args.db_type or 'nucl', parse_seqids=args.parse_seqids, verbose=args.verbose)

    blast_files = []
    if args.blast:
        blast_files = blaster.blast_run(
            input_path=args.query,
            db=args.name,
            db_dir=db_dir,
            evalue=float(args.evalue) or 1e-49,
            num_threads=int(args.num_threads),
            outfmt=args.outfmt or '6',
            output_dir=args.output_dir,
            verbose=args.verbose)

    go_files = []
    if args.go_slim_list:
        go_files = blaster.go_mart_to_go_slim_lists(
            go_file=args.go_slim_list, output_dir=args.output_dir)

    if args.split:
        blaster.split_blast_res_by_gos(
            blast_file_path=blast_files or args.in_blast_file,
            go_file_path=go_files or args.go,
            seqs_file_path=args.split,
            blast_outfmt=args.outfmt or 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            output_dir=args.output_dir,
            verbose=args.verbose)


def perform_stats(args: argparse.ArgumentParser) -> None:
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

    pcoa_path = ''
    if args.pcoa:
        pcoa_path = stats.pcoa(data_path=data_path, output_dir=args.output_dir, n_dim=args.n_dimensions,
                               plot=args.draw_plot, triangle=args.not_triangle, verbose=args.verbose)
    if pcoa_path:
        pcoa_path = pathlib.Path(pcoa_path)

    if args.anova:
        stats.anova(data_path=pcoa_path or data_path, groups_file=args.groups_file,
                    mode=args.mode, output_dir=args.output_dir, pcs=args.pc_number, ss_type=args.ss_type, triangle=args.not_triangle, verbose=args.verbose)

    if args.manova:
        stats.manova(data_path=pcoa_path or data_path,
                     groups_file=args.groups_file, mode=args.mode, output_dir=args.output_dir, pcs=args.pc_number, verbose=args.verbose)


def perform_mash(args: argparse.ArgumentParser) -> None:
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
        '-o', '--output_dir', help='output directory for the results')

    blaster_parser.add_argument(
        '-outfmt', '--outfmt', help='output format for BLAST results. Default: 6. Uses BLAST+ format codes.')

    blaster_parser.add_argument(
        '-v', '--verbose', action='store_true', help='verbose output')

    blaster_parser.add_argument(
        '-n', '--name', help='name of the BLAST database to create or use')

    # BLASTER - CREATING BLAST DATABASE
    blaster_parser.add_argument(
        '-cdb', '--create_db', action='store_true', help='create BLAST database from the -db_fasta FASTA file')

    blaster_parser.add_argument(
        '-dbt', '--db_type', help='type of the BLAST database to create (nucl or prot). Default: nucl')

    blaster_parser.add_argument('--parse_seqids', action='store_false',
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

    blaster_parser.add_argument('-e', '--evalue', type=float,
                                help='e-value threshold for BLAST search. Default: 10e-50')

    blaster_parser.add_argument('--num_threads', type=int, default=4,
                                help='number of threads to use for BLAST search. Default: 4')

    # CREATING GO SLIM LISTS FROM GO MART FILE
    blaster_parser.add_argument(
        '-gsl', '--go_slim_list', help='create GO slim lists from provided GO Mart file')

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
    stats_parser.add_argument(
        '-d', '--draw_plot', nargs=2, help='draw PCoA plot for chosen PCs. Two intigers required.')
    stats_parser.add_argument(
        '-g', '--groups_file', help='location of the file containing information on grouping for ANOVA or MANOVA. Required if -a or -ma was selected')
    stats_parser.add_argument('-ma', '--manova', action='store_true',
                              help='perfom MANOVA analysis on selected files')
    stats_parser.add_argument(
        '-mo', '--mode', default='n', help='select mode of ANOVA to perform. Should be either \'n\' (to perform ANOVA on all parameters), an integer (for m-way ANOVA where first m columns from the groups_file will be selected) or \'repeat\' for ANOVA with repeats. Defaults to \'n\'')
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
        if (bool(args.manova) or bool(args.anova)) ^ bool(args.groups_file):
            global_parser.error(
                '--anova or --manova and --groups_file must be given together!')

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
    # TODO maybe implement own parser for input from file and pass the correctly parsed args from stdout to the main function?
    main()
