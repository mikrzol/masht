import subprocess
import argparse
import pathlib

import mash_copy as mash


def main():
    global_parser = argparse.ArgumentParser(prog='masht',
                                            description='MASH distance toolkit',
                                            fromfile_prefix_chars='@')

    subparsers = global_parser.add_subparsers(
        title='subcommands', help='available subcommands')

    # mash subcommand
    mash_parser = subparsers.add_parser('mash', help='use the mash module')
    basic = mash_parser.add_argument_group('basic use')
    basic.add_argument('in_d',
                       help='location of 1) the folder with FASTQ or FASTA \
                           files or 2) the file with names of selected files \
                               (names are relative to package location)')
    detailed = mash_parser.add_argument_group('detailed use')
    detailed.add_argument('-d', '--distance', action='store_true',
                          help='calculate distance between selected files')
    detailed.add_argument('-s', '--sketch', action='store_true',
                          help='generate mash sketches of selected files')
    detailed.add_argument('-o', '--output_dir', default='./',
                          help='location of the output directory (default: ".")')
    detailed.add_argument('-v', '--verbose', action='store_true',
                          help='add more descriptions of performed actions')
    detailed.add_argument(
        '-m', '--mash', help='run mash with specified params', type=str)

    # parse args
    args = global_parser.parse_args()

    # variables
    data_path = pathlib.Path(args.in_d)

    # TESTING - remove paths later
    # testing_paths = ['../server_data/', '../mash-Linux64-v2.3/']
    bin_paths = ['bin/']

    # check if output dir exists and create one if necessary
    args.output_dir = args.output_dir.strip()
    pathlib.Path(f'{args.output_dir}').mkdir(
        parents=True, exist_ok=True)

    # mash_copy run()

    mash.run(vars(args))

    '''
    # dist
    if args.distance:
        mash.dist(bin_paths=bin_paths, data_path=data_path,
                  output_path=args.output_dir, verbose=args.verbose)

    # sketch
    if args.sketch:
        mash.sketch(bin_paths=bin_paths,
                    data_path=data_path, output_path=args.output_dir,
                    verbose=args.verbose)

    # mash (fully custom)
    if args.mash:
        subprocess.run(args.mash.split())
    '''


if __name__ == '__main__':
    main()
