import pathlib
from mash import _get_files
import pandas as pd


def _get_full_dist_matrix(df: pd.DataFrame) -> pd.DataFrame:
    from numpy import unique
    seqs = unique([df['seq_A']] + [df['seq_B']])
    actual_df = pd.DataFrame(index=seqs, columns=seqs)

    for _, row in df.iterrows():
        actual_df.loc[row['seq_A'], row['seq_B']] = row['distance']
        actual_df.loc[row['seq_B'], row['seq_A']] = row['distance']
        actual_df.loc[row['seq_A'], row['seq_A']] = 0
        actual_df.loc[row['seq_B'], row['seq_B']] = 0

    return actual_df


def pcoa(data_path: pathlib.Path, output_dir: pathlib.Path, verbose: bool = False) -> None:
    import pandas as pd
    from skbio.stats.ordination import pcoa as skb_pcoa

    files = _get_files(data_path)
    for file in files:
        df = pd.read_csv(file, sep='\t')
        df = _get_full_dist_matrix(df)

        res = skb_pcoa(df, number_of_dimensions=2)

        if verbose:
            print(f'********** {file.name.split(".")[0]} **********')
            print('========== Coordinates of samples in the ordination space: ==========')
            print(res.samples)
            print('\n========== Eigenvalues: ==========')
            print(res.eigvals)
            print('\n========== Proportion explained: ==========')
            print(res.proportion_explained)

        # create results file
        with open(f'{output_dir}/{file.name.split(".")[0]}_pcoa_results.txt', 'w') as out_fh:
            out_fh.write(
                '========== Coordinates of samples in the ordination space: ==========\n')
            out_fh.write(res.samples.to_string())
            out_fh.write('\n\n========== Eigenvalues: ==========\n')
            out_fh.write(res.eigvals.to_string())
            out_fh.write('\n\n========== Proportion explained: ==========\n')
            out_fh.write(res.proportion_explained.to_string())


pcoa(pathlib.Path('test_outputs/old_sketches_triangle.tsv'),
     pathlib.Path('test_outputs/'), True)
