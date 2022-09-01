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

combined_filtered = combined.loc[(combined['ORF_ID'].duplicated(keep=False) == False) | (combined['reference_db'] == 'BLDB'), :]

sorted = combined_filtered.sort_values(by=['Sample'])

sorted.to_csv(args.output, sep="\t", index=None)
