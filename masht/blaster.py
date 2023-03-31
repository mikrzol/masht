import pandas as pd
import pathlib


def query_for_go_terms() -> None:
    # TODO implement querying for latest GO plant terms file from GO server
    pass


def go_mart_to_go_slim_lists(go_file: str, output_dir: str) -> list[str]:
    # create output dir
    pathlib.Path(f'{output_dir}/go_lists').mkdir(
        parents=True, exist_ok=True)

    # read in go file
    go_df = pd.read_csv(go_file)

    # group by go slim terms (assuming last column is go slim term)
    grouped = go_df.groupby(go_df.columns[-1])

    # save to appropriate csv files
    files = []
    for name, group in grouped:
        group.to_csv(f'{output_dir}/go_lists/{name}.csv', index=False)
        files.append(f'{output_dir}/go_lists/{name}.csv')

    return files
