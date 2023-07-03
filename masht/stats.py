import pathlib
from mash import _get_files
import pandas as pd


def analyze_all(data_path: pathlib.Path, mode: str, groups_file: str, output_dir: pathlib.Path, formula: str or None, anova_manova_mode: str = 'n', pcs: int = 4, verbose: bool = False, plot: list[str] = ['1', '2'], ss_type: int = 2, triangle: bool = True, n_dim: int or None = None) -> None:
    """Perform standardized full analysis of all files in the data_path directory

    Args:
        data_path (pathlib.Path): location of directory with subdirs with triangle (.tsv) files from mash.analyze_all
        mode (str): what type of analysis to perform. Can be either 'anova' or 'manova'
        groups_file (str): path to the file with group descriptions. The first column of this file should contain names of observations.
        output_dir (pathlib.Path): where to save the results
        anova_manova_mode (str, optional): How many parameters from the groups_file to consider. Defaults to 'n', which means to consider all of them.
        pcs (int, optional): How many PCs to consider when performing ANOVA/MANOVA. Defaults to 4.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.
        plot (list[str], optional): which PCs to plot when performing PCoA. Defaults to ['1', '2']. If not provided, no plots will be generated.
        ss_type (int, optional): which sum of squares type to use in ANOVA analysis. Defaults to 2. Look up statsmodels.stats.anova.anova_lm for more details.
        triangle (bool, optional): whether the input file has the triangle form. Defaults to True.
        n_dim (int or None, optional): number of dimensions to use in PCoA. Defaults to None (meaning equal to number of observations).
    """

    from joblib import Parallel, delayed
    from collections import Counter

    def _mp_analyze(file: pathlib.Path):
        if verbose:
            print(f'Performing PCoA on {str(file.parent).split("/")[-1]}...')
        try:
            pcoa_path = pathlib.Path(pcoa(data_path=file, output_dir=file.parent,
                                          n_dim=n_dim, plot=plot, triangle=triangle, verbose=False))
        except Exception as e:
            print(f"Error occurred for {str(file)}: {str(e)}")
            return

        if mode == 'anova':
            if verbose:
                print(
                    f'Performing ANOVA on {str(file.parent).split("/")[-1]}...')
            try:
                anova(data_path=pcoa_path, groups_file=groups_file, output_dir=file.parent, formula=formula,
                      anova_manova_mode=anova_manova_mode, pcs=pcs, ss_type=ss_type, triangle=triangle, verbose=False)
            except Exception as e:
                print(f"Error occurred for {str(file)}: {str(e)}")
                return

        else:
            if verbose:
                print(
                    f'Performing MANOVA on {str(file.parent).split("/")[-1]}...')
            try:
                manova(data_path=pcoa_path, groups_file=groups_file, formula=formula,
                       output_dir=file.parent, anova_manova_mode=anova_manova_mode, pcs=pcs, verbose=False)
            except Exception as e:
                print(f"Error occurred for {str(file)}: {str(e)}")
                return

    # get number of filtered_* files in each subdir to determine which subdirs to analyze
    lst = [x.parent for x in pathlib.Path(
        data_path).rglob('*filtered_*.fasta')]

    counted_lst = Counter(lst)

    separator = ",\n\t"
    print(
        f'masht stats analyze_all: Omitting {sum([1 for x in counted_lst.values() if x < max(counted_lst.values())])} subdirectories:\n\t{separator.join([str(k) for k, v in counted_lst.items() if v < max(counted_lst.values())])}\nwith less than {max(counted_lst.values())} *filtered_*.fasta files.\n')

    my_dict = {k: v for k, v in counted_lst.items() if v ==
               max(counted_lst.values())}

    subdirs = [x / 'sketches_triangle.tsv' for x in my_dict.keys()]

    Parallel(n_jobs=-1)(delayed(_mp_analyze)(subdir)
                        for subdir in subdirs)


def manova(data_path: pathlib.Path, groups_file: str, output_dir: pathlib.Path, formula: str or None, anova_manova_mode: str = 'n', pcs: int = 4, verbose: bool = False) -> None:
    """perform MANOVA analyses of selected files

    Args:
        data_path (pathlib.Path): location of the input file (with PCoA coordinates). Automatically inferred if pcoa step was performed in the same command. 
        groups_file (str): location of the file with group descriptions. The first column of this file should contain names 
        output_dir (pathlib.Path): output directory location
        anova_manova_mode (str, optional): How many parameters from the groups_file to consider. Defaults to 'n', which means to consider all of them.
        pcs (int, optional): how many PCs to consider. Defaults to 4.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.
    """

    # TODO add checking for assumptions correctness with
    # 1. generalized shapiro-wilk's for normality
    # 2. box's m test for homogeneity of covariances
    # 3. checking for multicollinearity
    # 4. checking for linear relationships of the dependent variables

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
        # remove prefix from index - for full analysis (analyze_all)
        df.index = df.index.str.removeprefix(
            'filtered_').str.split('.').str.get(0)

        groups = pd.read_csv(groups_file, sep='\t', index_col=0)

        # select columns with non-zero values only
        df = df[df.columns[~(df == 0).all()]]

        # trim to the num of columns that are non zero
        if pcs > len(df.columns):
            pcs = len(df.columns)

        manova = pd.merge(groups, df, left_index=True, right_index=True)

        # get number of params to consider
        if anova_manova_mode.isnumeric():
            anova_manova_mode = int(anova_manova_mode)
        else:
            anova_manova_mode = len(groups.columns)

        # convert Pandas df to R df
        with (ro.default_converter + pandas2ri.converter).context():
            m_df = ro.conversion.get_conversion().py2rpy(manova)

        # create formula and create model
        from rpy2.robjects import Formula
        if not formula:
            form = " * ".join(groups.columns[:anova_manova_mode])
        else:
            form = formula
        man_formula = Formula(
            f'cbind({",".join(df.columns[0:pcs])}) ~ {form}')
        man = stats.manova(formula=man_formula, data=m_df)

        # convert results to Pandas df
        df_list = []
        with (ro.default_converter + pandas2ri.converter).context():
            for test in ['Wilks', 'Pillai', 'Hotelling-Lawley', 'Roy']:
                temp = ro.conversion.get_conversion().rpy2py(
                    broom.tidy_manova(man, test=test))
                temp.rename({temp.columns[2]: 'value'}, axis=1, inplace=True)
                temp['test'] = test
                temp.set_index(['test', 'term'], inplace=True)
                df_list.append(temp)
        big_df = pd.concat(df_list)

        # save results to file
        if not formula:
            big_df.to_csv(
                f'{output_dir}/manova_{pcs}_PCs_{anova_manova_mode}_params.tsv', sep='\t')
        else:
            big_df.to_csv(
                f'{output_dir}/manova_{pcs}_PCs_{formula.replace(" ", "")}.tsv', sep='\t')

        # report results if verbose
        if verbose:
            print('========= MANOVA results: =========')
            print(big_df)

    ''' implementation with Python's statsmodels
        from statsmodels.multivariate.manova import MANOVA

        # create formula and create model
        formula = f'{" + ".join(df.columns[0:pcs])} ~ {" * ".join(groups.columns[:anova_manova_mode])}'
        fit = MANOVA.from_formula(formula, data=manova)

        # report results and save them to file
        if verbose:
            print()
            print(fit.mv_test().summary_frame)

        fit.mv_test().summary_frame.to_csv(
            f'{output_dir}/manova_{pcs}_PCs_{anova_manova_mode}_params.tsv', sep='\t')
    '''


def anova(data_path: pathlib.Path, groups_file: str, output_dir: pathlib.Path, formula: str or None, anova_manova_mode: str = 'n', pcs: int = 4, ss_type: int = 2, triangle: bool = True, verbose: bool = False) -> None:
    """perform ANOVA (either with repeats or n-way) on data in data_path file

    Args:
        data_path (pathlib.Path): location of the input file (with PCoA coordinates). Automatically inferred if pcoa step was performed in the same command.
        groups_file (str): location of the file with group descriptions. The first column of this file should contain names 
        output_dir (pathlib.Path): output directory location
        anova_manova_mode (str, optional): which type of ANOVA to perform. Use 'repeat' for repeated ANOVA, int to perform n-way ANOVA or leave as default to perform ANOVA on all parameters in the groups_file. Defaults to 'n' (for all parameters).
        pcs (int, optional): number of PCs to perform ANOVA analyses on. Defaults to 4.
        ss_type (int, optional): the type of sum of squares to use (ANOVA-specific). Implemented in statsmodels.anova_lm. Defaults to 2.
        triangle (bool, optional): (not used) whether input is in the form of a triangle. Defaults to True.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.
    """

    files = _get_files(data_path)

    for file in files:
        # print('working on file:', file)
        df = pd.read_csv(file, sep='\t', index_col=0)
        # remove prefix from index - for full analysis (analyze_all)
        df.index = df.index.str.removeprefix(
            'filtered_').str.split('.').str.get(0)

        # select columns with non-zero values only
        df = df[df.columns[~(df == 0).all()]]

        # trim to the num of columns that are non zero
        if pcs > len(df.columns):
            pcs = len(df.columns)

        groups = pd.read_csv(groups_file, sep='\t', index_col=0)

        if anova_manova_mode == 'repeat':
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
            n = int(anova_manova_mode) if anova_manova_mode.isnumeric() else len(
                groups.columns)

            # select columns
            groups_selected = groups.iloc[:, :n]

            aov_df = pd.merge(groups_selected, df,
                              left_index=True, right_index=True)

            # perform n-way ANOVAs
            # values in brackets get PC1-PCn (parameters are before PCs in aov_df)
            for col in aov_df.columns[len(groups_selected.columns):len(groups_selected.columns)+pcs]:
                # perform full ANOVA on selected factors (all interactions included)

                if not formula:
                    formula = f'{" * ".join(groups_selected.columns)}'

                model = ols(f'{col} ~ {formula}', data=aov_df).fit()

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


def plot_pcoa(res, names: list[str], output_dir: pathlib.Path, chosen_pcs: list[str] = ['1', '2']) -> None:
    """plot skbio.(...).pcoa results and save it to file

    Args:
        res (_type_): OrdinationResults object created by skb_pcoa function
        names (list[str]): list of names of observations
        output_dir (pathlib.Path): output location
        chosen_pcs (list[str], optional): which PCs to plot. Defaults to ['1', '2'].
    """

    # matplotlib implementation
    import matplotlib.pyplot as plt
    chosen_pcs = ['PC' + pc for pc in chosen_pcs]

    plt.figure(figsize=(10, 8), dpi=200)
    plt.plot(res.samples[chosen_pcs[0]],
             res.samples[chosen_pcs[1]], 'o')
    plt.grid(color='lightgrey')
    for name, (i, pc) in zip(names, res.samples[chosen_pcs].iterrows()):
        plt.annotate(name, pc, xytext=(10, -5),
                     textcoords='offset points', color='darkslategrey', annotation_clip=True)
    plt.title('PCoA ordination')
    plt.xlabel(
        f'{chosen_pcs[0]} ({round(res.proportion_explained[chosen_pcs[0]]*100,2)}% variance explained)')
    plt.ylabel(
        f'{chosen_pcs[1]} ({round(res.proportion_explained[chosen_pcs[1]]*100,2)}% variance explained)')
    plt.savefig(f'{str(output_dir)}/pcoa_plot.png', bbox_inches='tight')

    ''' # plot using plotnine
    from plotnine import ggplot, aes, geom_point, geom_text, labs, scale_x_continuous

    names = res.samples.index
    chosen_pcs = ['PC' + pc for pc in chosen_pcs]
    nudge_val = res.samples['PC1'].max() * 0.03

    plot = (
        ggplot(res.samples) +
        aes(x=chosen_pcs[0], y=chosen_pcs[1]) +
        geom_text(aes(label=names), nudge_x=nudge_val, nudge_y=nudge_val) +
        geom_point() +
        labs(title=f'PCoA ordination of {" and ".join(chosen_pcs)}',
             x=f'{chosen_pcs[0]} ({round(res.proportion_explained[chosen_pcs[0]]*100,2)}% variance explained)',
             y=f'{chosen_pcs[1]} ({round(res.proportion_explained[chosen_pcs[1]]*100,2)}% variance explained)')
        + scale_x_continuous(expand=(0.1, 0))
    )

    plot.save(f'{str(output_dir)}/pcoa_plot.png',
              width=10, height=10, dpi=400, verbose=False)
    '''


def pcoa(data_path: pathlib.Path, output_dir: pathlib.Path, n_dim: int or None = None, plot: list[str] = [], triangle: bool = True, verbose: bool = False) -> str:
    """perform PCoA of data obtained in the mash.triangle function

    Args:
        data_path (pathlib.Path): location of input file
        output_dir (pathlib.Path): output location
        n_dim (int or None, optional): number of dimensions to use in PCoA. Defaults to None (meaning equal to number of observations).
        plot (list[str], optional): which PCs to plot if any. Defaults to [] (meaning: don't plot).
        triangle (bool, optional): whether the input file has the triangle form. Defaults to True.
        verbose (bool, optional): whether to increase verbosity. Defaults to False.


    Returns:
        str: location of the last generated pcoa_coords.tsv file. Used to automatically perform the anova step on such file if performing ANOVA was specified in the command. 
    """
    import pandas as pd
    from skbio.stats.ordination import pcoa as skb_pcoa
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

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
                    res.samples.drop(
                        res.samples.columns[-1], axis=1, inplace=True)
                    res.eigvals = res.eigvals[:-1]
                    res.proportion_explained = res.proportion_explained[:-1]
            res.eigvals.rename('eigenval', inplace=True)
            res.proportion_explained.rename(
                'proportion_explained', inplace=True)

            if verbose:
                print()
                print(f'********** {file.name.split(".")[0]} **********')
                print(
                    '========== Coordinates of samples in the ordination space: ==========')
                print(res.samples)
                print('\n========== Eigenvalues: ==========')
                print(res.eigvals)
                print('\n========== Proportion explained: ==========')
                print(res.proportion_explained)

            # plot
            if plot:
                if triangle:
                    plot_pcoa(res=res, names=names,
                              output_dir=output_dir, chosen_pcs=plot)
                else:
                    plot_pcoa(res=res, names=df.columns,
                              output_dir=output_dir, chosen_pcs=plot)

            # create results file
            res.samples.to_csv(
                f'{output_dir}/{file.name.split(".")[0]}_pcoa_coords.tsv', sep='\t')
            last_result = f'{output_dir}/{file.name.split(".")[0]}_pcoa_coords.tsv'
            res.eigvals.to_csv(
                f'{output_dir}/{file.name.split(".")[0]}_pcoa_eigenvals.tsv', sep='\t')
            res.proportion_explained.to_csv(
                f'{output_dir}/{file.name.split(".")[0]}_pcoa_proportions.tsv', sep='\t')

        return last_result
