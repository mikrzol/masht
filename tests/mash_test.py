from masht.mash import sketch
import pathlib

# variables
data_dir = pathlib.Path('../server_data/seqs')
data_file = pathlib.Path('tests/files.txt')
testing_paths = ['../server_data/', '../mash-Linux64-v2.3/']

# basic sketch
sketch(paths=testing_paths, data_path=data_file,
       output_path='test_outputs')
