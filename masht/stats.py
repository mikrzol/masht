import pathlib
from mash import _get_files
import pandas as pd


def _get_full_dist_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """convert triangle matrix to full distance matrix

    Args:
        df (pd.DataFrame): df created by reading in results from mash.triangle function 

    Returns:
        pd.DataFrame: converted DataFrame
    """
    from numpy import unique
    seqs = unique([df['seq_A']] + [df['seq_B']])
    actual_df = pd.DataFrame(index=seqs, columns=seqs)

    for _, row in df.iterrows():
        actual_df.loc[row['seq_A'], row['seq_B']] = row['distance']
        actual_df.loc[row['seq_B'], row['seq_A']] = row['distance']
        actual_df.loc[row['seq_A'], row['seq_A']] = 0
        actual_df.loc[row['seq_B'], row['seq_B']] = 0

    return actual_df


def plot_pcoa(res, names: list[str], output_dir: pathlib.Path) -> None:
    """plot skbio.(...).pcoa results and save it to file

    Args:
        res (_type_): OrdinationResults object created by skb_pcoa function
        names (list[str]): list of names of observations
        output_dir (pathlib.Path): output location
    """
    import matplotlib.pyplot as plt

    plt.plot(res.samples['PC1'], res.samples['PC2'], 'o')
    plt.grid(color='lightgrey')
    for name, (i, pc) in zip(names, res.samples[['PC1', 'PC2']].iterrows()):
        plt.annotate(name, pc, xytext=(10, -5),
                     textcoords='offset points', color='darkslategrey', annotation_clip=True)
    plt.title('PCoA ordination')
    plt.xlabel(
        f'PC1 ({round(res.proportion_explained["PC1"]*100,2)}% variance explained)')
    plt.ylabel(
        f'PC2 ({round(res.proportion_explained["PC2"]*100,2)}% variance explained)')
    plt.savefig(f'{str(output_dir)}/pcoa_plot.png', bbox_inches='tight')


def pcoa(data_path: pathlib.Path, output_dir: pathlib.Path, n_dim: int or None = None, plot: bool = False, verbose: bool = False) -> None:
    """perform PCoA of data obtained in the mash.triangle function

    Args:
        data_path (pathlib.Path): location of input triangle file
        output_dir (pathlib.Path): output location
        n_dim (intorNone, optional): number of dimensions to use in PCoA. Defaults to None (meaning equal to number of observations).
        verbose (bool, optional): whether to increase verbosity. Defaults to False.
    """
    import pandas as pd
    from skbio.stats.ordination import pcoa as skb_pcoa

    files = _get_files(data_path)
    for file in files:
        df = pd.read_csv(file, sep='\t')
        df = _get_full_dist_matrix(df)

        res = skb_pcoa(df, number_of_dimensions=n_dim or len(df))

        if verbose:
            print(f'********** {file.name.split(".")[0]} **********')
            print('========== Coordinates of samples in the ordination space: ==========')
            print(res.samples)
            print('\n========== Eigenvalues: ==========')
            print(res.eigvals)
            print('\n========== Proportion explained: ==========')
            print(res.proportion_explained)

        # plot
        if plot:
            names = [x.split('/')[-1] for x in df.columns]
            plot_pcoa(res=res, names=names, output_dir=output_dir)

        # create results file
        res.samples.to_csv(
            f'{output_dir}/{file.name.split(".")[0]}_pcoa_coords.csv')
        res.eigvals.to_csv(
            f'{output_dir}/{file.name.split(".")[0]}_pcoa_eigenvals.csv')
        res.proportion_explained.to_csv(
            f'{output_dir}/{file.name.split(".")[0]}_pcoa_proportions.csv')
