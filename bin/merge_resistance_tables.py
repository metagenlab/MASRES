#!/opt/conda/envs/python-r/bin/python
import argparse
import os
from typing import NamedTuple
import pandas 

class Args(NamedTuple):
    """ Command-line arguments """
    tsv_files: str
    output: str

def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Merge amr predictions for all samples in a run',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('tsv_files',
                        metavar='tsv files',
                        help='tsv files with AMR predictions',
			nargs='+')
    parser.add_argument('-o',
                        '--out',
                        metavar='output file',
                        type=str,
                        help='name of the output file',
                        default='merged_table_resistance.tsv',
                        )

    args = parser.parse_args()

    return Args(args.tsv_files, args.out)

args = get_args()

df_list = [pandas.read_csv(table, sep="\t", index_col=None, header=0) for table in args.tsv_files]

combined = pandas.concat(df_list)

list_of_df = [v for k, v in combined.groupby("ORF_ID")] #split df by ORF ID
for i, df in enumerate(list_of_df):
    if df.shape[0] > 1: #If predictions were made by both BLDB and RGI for the same ORF
        list_of_df[i] = df[df['Percent_coverage']==df['Percent_coverage'].max()] #Keep the one with the highest perc coverage
final_df = pandas.concat(list_of_df)

final_df["Sample"] = final_df["Sample"].astype(str)

sorted = final_df.sort_values(by=['Sample'])

sorted.to_csv(args.output, sep="\t", index=None)
