import pathlib
from mash import _get_files
import pandas as pd


def manova(data_path: pathlib.Path, groups_file: str, output_dir: pathlib.Path, mode: str = 'n', pcs: int = 4, verbose: bool = False) -> None:
    """perform MANOVA analyses of selected files

    Args:
        data_path (pathlib.Path): location of the input file (with PCoA coordinates). Automatically inferred if pcoa step was performed in the same command. 
        groups_file (str): location of the file with group descriptions. The first column of this file should contain names 
        output_dir (pathlib.Path): output directory location
        mode (str, optional): How many parameters from the groups_file to consider. Defaults to 'n', which means to consider all of them.
        pcs (int, optional): how many PCs to consider. Defaults to 4.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.
    """

    # using rpy2 to perform MANOVA because Python's statsmodels gives incorrect results

    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    import rpy2.robjects as ro
    stats = importr('stats')
    broom = importr('broom')

    # supress meaningless errors
    import rpy2.rinterface_lib.callbacks
    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = lambda *args: None

    files = _get_files(data_path)

    for file in files:
        df = pd.read_csv(file, sep='\t', index_col=0)
        groups = pd.read_csv(groups_file, sep='\t', index_col=0)

        # select columns with non-zero values only
        df = df[df.columns[~(df == 0).all()]]

        manova = pd.merge(groups, df, left_index=True, right_index=True)

        # get number of params to consider
        if mode.isnumeric():
            mode = int(mode)
        else:
            mode = len(groups.columns)

        # convert Pandas df to R df
        with (ro.default_converter + pandas2ri.converter).context():
            m_df = ro.conversion.get_conversion().py2rpy(manova)

        # create formula and create model
        from rpy2.robjects import Formula
        formula = Formula(
            f'cbind({",".join(df.columns[0:pcs])}) ~ {" * ".join(groups.columns[:mode])}')
        man = stats.manova(formula=formula, data=m_df)

        # convert results to Pandas df
        df_list = []
        with (ro.default_converter + pandas2ri.converter).context():
            for test in ['Wilks', 'Pillai', 'Hotelling-Lawley', 'Roy']:
                temp = ro.conversion.get_conversion().rpy2py(broom.tidy_manova(man, test=test))
                temp.rename({temp.columns[2]: 'value'}, axis=1, inplace=True)
                temp['test'] = test
                temp.set_index(['test', 'term'], inplace=True)
                df_list.append(temp)
        big_df = pd.concat(df_list)

        # save results to file
        big_df.to_csv(
            f'{output_dir}/manova_{pcs}_PCs_{mode}_params.tsv', sep='\t')

        # report results if verbose
        if verbose:
            print('========= MANOVA results: =========')
            print(big_df)

    ''' implementation with Python's statsmodels
        from statsmodels.multivariate.manova import MANOVA

        # create formula and create model
        formula = f'{" + ".join(df.columns[0:pcs])} ~ {" * ".join(groups.columns[:mode])}'
        fit = MANOVA.from_formula(formula, data=manova)

        # report results and save them to file
        if verbose:
            print()
            print(fit.mv_test().summary_frame)

        fit.mv_test().summary_frame.to_csv(
            f'{output_dir}/manova_{pcs}_PCs_{mode}_params.tsv', sep='\t')
    '''


def anova(data_path: pathlib.Path, groups_file: str, output_dir: pathlib.Path, mode: str = 'n', pcs: int = 4, ss_type: int = 2, triangle: bool = True, verbose: bool = False) -> None:
    """perform ANOVA (either with repeats or n-way) on data in data_path file

    Args:
        data_path (pathlib.Path): location of the input file (with PCoA coordinates). Automatically inferred if pcoa step was performed in the same command.
        groups_file (str): location of the file with group descriptions. The first column of this file should contain names 
        output_dir (pathlib.Path): output directory location
        mode (str, optional): which type of ANOVA to perform. Use 'repeat' for repeated ANOVA, int to perform n-way ANOVA or leave as default to perform ANOVA on all parameters in the groups_file. Defaults to 'n' (for all parameters).
        pcs (int, optional): number of PCs to perform ANOVA analyses on. Defaults to 4.
        ss_type (int, optional): the type of sum of squares to use (ANOVA-specific). Implemented in statsmodels.anova_lm. Defaults to 2.
        triangle (bool, optional): (not used) whether input is in the form of a triangle. Defaults to True.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.
    """
    files = _get_files(data_path)

    for file in files:
        df = pd.read_csv(file, sep='\t', index_col=0)
        groups = pd.read_csv(groups_file, sep='\t', index_col=0)

        if mode == 'repeat':
            from statsmodels.stats.anova import AnovaRM

            # get dataframe
            aov_df = pd.merge(groups, df, left_index=True, right_index=True)

            # prepare results container DF
            aovs = pd.DataFrame(
                columns=['F Value', 'Num DF', 'Den DF', 'Pr > F'])

            # perform repeated ANOVA
            group_names = aov_df.columns[:2]
            for col in aov_df.columns[2:]:
                temp = AnovaRM(data=aov_df, depvar=col,
                               subject=group_names[0], within=group_names[1]).fit()
                temp.anova_table.index = [col]
                aovs = pd.concat([aovs, temp.anova_table])

            if verbose:
                print()
                print('ANOVA repeated:')
                print(aovs)

            aovs.to_csv(f'{output_dir}/aov_repeated.tsv', sep='\t')

        else:
            import statsmodels.api as sm
            from statsmodels.formula.api import ols

            # get number of parameters to perform ANOVA on
            n = int(mode) if mode.isnumeric() else len(groups.columns)

            # select columns
            groups_selected = groups.iloc[:, :n]

            aov_df = pd.merge(groups_selected, df,
                              left_index=True, right_index=True)

            # perform n-way ANOVAs
            # values in brackets get PC1-PCn (parameters are before PCs in aov_df)
            for col in aov_df.columns[len(groups_selected.columns):len(groups_selected.columns)+pcs]:
                # perform full ANOVA on selected factors (all interactions included)
                model = ols(
                    f'{col} ~ {" * ".join(groups_selected.columns)}', data=aov_df).fit()

                if verbose:
                    print()
                    print(f'ANOVA {len(groups_selected.columns)}-way ({col}):')
                    print(sm.stats.anova_lm(model, typ=ss_type))
                    print()

                sm.stats.anova_lm(model, typ=ss_type).to_csv(
                    f'{output_dir}/aov_{n}-way_{col}.tsv', sep='\t')


def _get_full_dist_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """convert triangle matrix to full distance matrix

    Args:
        df (pd.DataFrame): df created by reading in results from mash.triangle function 

    Returns:
        pd.DataFrame: converted DataFrame
    """
    from numpy import unique
    col_1, col_2, col_3 = df.columns
    seqs = unique([df[col_1]] + [df[col_2]])
    actual_df = pd.DataFrame(index=seqs, columns=seqs)

    for _, row in df.iterrows():
        actual_df.loc[row[col_1], row[col_2]] = row[col_3]
        actual_df.loc[row[col_2], row[col_1]] = row[col_3]
        actual_df.loc[row[col_1], row[col_1]] = 0
        actual_df.loc[row[col_2], row[col_2]] = 0

    return actual_df


def plot_pcoa(res, names: list[str], output_dir: pathlib.Path) -> None:
    """plot skbio.(...).pcoa results and save it to file

    Args:
        res (_type_): OrdinationResults object created by skb_pcoa function
        names (list[str]): list of names of observations
        output_dir (pathlib.Path): output location
    """
    import matplotlib.pyplot as plt

    # TODO implement option to specify which PCs to plot

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


def pcoa(data_path: pathlib.Path, output_dir: pathlib.Path, n_dim: int or None = None, plot: bool = False, triangle: bool = True, verbose: bool = False) -> str:
    """perform PCoA of data obtained in the mash.triangle function

    Args:
        data_path (pathlib.Path): location of input file
        output_dir (pathlib.Path): output location
        n_dim (intorNone, optional): number of dimensions to use in PCoA. Defaults to None (meaning equal to number of observations).
        plot (bool, optional): whether to create a plot of results. Defaults to False.
        triangle (bool, optional): whether the input file has the triangle form. Defaults to True.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.

    Returns:
        str: location of the last generated pcoa_coords.tsv file. Used to automatically perform the anova step on such file if performing ANOVA was specified in the command. 
    """
    import pandas as pd
    from skbio.stats.ordination import pcoa as skb_pcoa

    last_result = ''  # location of file to perform ANOVA later if chosen
    files = _get_files(data_path)
    for file in files:
        if triangle:
            df = pd.read_csv(file, sep='\t')
            df = _get_full_dist_matrix(df)
        else:
            df = pd.read_csv(file, sep='\t', header=None)

        res = skb_pcoa(df, number_of_dimensions=n_dim or len(df))

        names = [x.split('/')[-1].split('.')[0] for x in df.columns]
        if triangle:
            # rename rows
            res.samples.index = names
            # clean results
            if not n_dim:
                res.samples.drop(res.samples.columns[-1], axis=1, inplace=True)
                res.eigvals = res.eigvals[:-1]
                res.proportion_explained = res.proportion_explained[:-1]
        res.eigvals.rename('eigenval', inplace=True)
        res.proportion_explained.rename(
            'proportion_explained', inplace=True)

        if verbose:
            print()
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
