import subprocess
import argparse
import pathlib

import mash
import stats


def perform_stats(args: argparse.ArgumentParser, bin_path: list[str], data_path: pathlib.Path) -> None:
    """perform the stats subcommand

    Args:
        args (argparse.ArgumentParser): args created in the main function
        bin_path (str): directory with the mash binary (obsolete)
        data_path (pathlib.Path): location of the dir/ path with files to work on
    """
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


def perform_mash(args: argparse.ArgumentParser, bin_path: str, data_path: pathlib.Path) -> None:
    """perform mash command

    Args:
        args (argparse.ArgumentParser): args created in the main function
        bin_path (str): directory with the mash binary
        data_path (pathlib.Path): location of the dir/ path with files to work on
    """
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

    # remove spaces from global_parser
    global_parser

    subparsers = global_parser.add_subparsers(
        title='subcommands', help='available subcommands', dest='subparser_name')

    # STATS SUBCOMMAND
    stats_parser = subparsers.add_parser('stats', help='use the stats module')
    stats_parser.add_argument('in_d',
                              help='location of 1) the folder with FASTQ or FASTA \
                           files or 2) the file with names of selected files \
                               (names are relative to package location)')
    stats_parser.add_argument(
        '-a', '--anova', action='store_true', help='perfom ANOVA on selected files')
    stats_parser.add_argument(
        '-d', '--draw_plot', action='store_true', help='draw plots for performed analyses')
    stats_parser.add_argument(
        '-g', '--groups_file', help='location of the file containing information on grouping for ANOVA. Required if -a was selected')
    stats_parser.add_argument('-ma', '--manova', action='store_true',
                              help='perfom MANOVA analysis on selected files')
    stats_parser.add_argument(
        '-mo', '--mode', default='n', help='select mode of ANOVA to perform. Should be either \'n\' (to perform ANOVA on all parameters), an integer (for m-way ANOVA where first m columns from the groups_file will be selected) or \'repeat\' for ANOVA with repeats.')
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
    args = global_parser.parse_args()
    if args.subparser_name == 'stats':
        if (bool(args.manova) or bool(args.anova)) ^ bool(args.groups_file):
            global_parser.error(
                '--anova or --manova and --groups_file must be given together!')

    # VARIABLES
    data_path = pathlib.Path(args.in_d)

    # TESTING - remove paths later
    bin_path = 'bin/'

    # check if output dir exists and create one if necessary
    args.output_dir = args.output_dir.strip()
    pathlib.Path(f'{args.output_dir}').mkdir(
        parents=True, exist_ok=True)

    # MASH
    # perform_mash(args=args, bin_path=bin_path, data_path=data_path)

    # STATS
    # perform_stats(args=args, bin_path=bin_path, data_path=data_path)

    # perform appropriate function
    args.func(args, bin_path, data_path)


if __name__ == '__main__':
    main()
