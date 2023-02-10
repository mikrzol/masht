import pathlib
from mash import _get_files
import pandas as pd


def anova(data_path: pathlib.Path, groups_file: str, mode: str, output_dir: pathlib.Path, verbose: bool = False):
    """perform ANOVA (either with repeats or two-way) on data in data_path file

    Args:
        data_path (pathlib.Path): location of file to perform ANOVA on
        groups_file (str): location of .tsv file with information on groups for ANOVA 
        mode (str): which type of ANOVA to perform
        output_dir (pathlib.Path): location of folder to put the results files into
        verbose (bool, optional): whether to report performed actions to console. Defaults to False.
    """
    files = _get_files(data_path)
    for file in files:
        # load in triangle file and clean it
        df = pd.read_csv(file, sep='\t', index_col=0)
        '''
        df.index = [x.split('/')[-1] for x in df.index]
        df.columns = [x.split('/')[-1] for x in df.columns]
        '''

        groups = pd.read_csv(groups_file, sep='\t')
        groups = groups.set_index('file')
        aov_df = pd.merge(groups, df, left_index=True, right_index=True)

        if mode == 'twoway':
            import statsmodels.api as sm
            from statsmodels.formula.api import ols

            # perform two-way ANOVAs
            for col in aov_df.columns[2:]:
                # TESTING - add '+ A:B' to anova formula below!!!
                model = ols(f'{col} ~ A + B', data=aov_df).fit()

                if verbose:
                    print(f'ANOVA two-way ({col}):')
                    print(sm.stats.anova_lm(model, typ=2))
                    print()

                sm.stats.anova_lm(model, typ=2).to_csv(
                    f'{output_dir}/aov_two-way_{col}.tsv', sep='\t')

        elif mode == 'repeat':
            from statsmodels.stats.anova import AnovaRM
            aovs = pd.DataFrame(
                columns=['F Value', 'Num DF', 'Den DF', 'Pr > F'])
            for col in aov_df.columns[2:]:
                temp = AnovaRM(data=aov_df, depvar=col,
                               subject='A', within=['B']).fit()
                temp.anova_table.index = [col]
                aovs = pd.concat([aovs, temp.anova_table])

            if verbose:
                print('ANOVA repeated:')
                print(aovs)

            aovs.to_csv(f'{output_dir}/aov_repeated.tsv', sep='\t')


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


def pcoa(data_path: pathlib.Path, output_dir: pathlib.Path, n_dim: int or None = None, plot: bool = False, triangle: bool = False, verbose: bool = False) -> str:
    """perform PCoA of data obtained in the mash.triangle function

    Args:
        data_path (pathlib.Path): location of input triangle file
        output_dir (pathlib.Path): output location
        n_dim (intorNone, optional): number of dimensions to use in PCoA. Defaults to None (meaning equal to number of observations).
        verbose (bool, optional): whether to increase verbosity. Defaults to False.
    """
    import pandas as pd
    from skbio.stats.ordination import pcoa as skb_pcoa

    last_result = ''  # location of file to perform ANOVA later if chosen
    files = _get_files(data_path)
    for file in files:
        df = pd.read_csv(file, sep='\t', header=None)
        if triangle:
            df = pd.read_csv(file, sep='\t')
            df = _get_full_dist_matrix(df)

        res = skb_pcoa(df, number_of_dimensions=n_dim or len(df))

        if triangle:
            # rename rows
            res.samples.index = [x.split('/')[-1] for x in df.columns]
            # clean results
            if not n_dim:
                res.samples.drop(res.samples.columns[-1], axis=1, inplace=True)
                res.eigvals = res.eigvals[:-1]
                res.proportion_explained = res.proportion_explained[:-1]
        res.eigvals.rename('eigenval', inplace=True)
        res.proportion_explained.rename(
            'proportion_explained', inplace=True)

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
            if triangle:
                names = [x.split('/')[-1] for x in df.columns]
                plot_pcoa(res=res, names=names, output_dir=output_dir)
            else:
                plot_pcoa(res=res, names=df.columns, output_dir=output_dir)

        # create results file
        res.samples.to_csv(
            f'{output_dir}/{file.name.split(".")[0]}_pcoa_coords.tsv', sep='\t')
        last_result = f'{output_dir}/{file.name.split(".")[0]}_pcoa_coords.tsv'
        res.eigvals.to_csv(
            f'{output_dir}/{file.name.split(".")[0]}_pcoa_eigenvals.tsv', sep='\t')
        res.proportion_explained.to_csv(
            f'{output_dir}/{file.name.split(".")[0]}_pcoa_proportions.tsv', sep='\t')

    return last_result
