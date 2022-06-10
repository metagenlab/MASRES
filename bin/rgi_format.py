#!/opt/conda/envs/python-r/bin/python
import argparse
from ast import Str
import os
from typing import NamedTuple
import pandas 
from Bio import SeqIO

class Args(NamedTuple):
    """ Command-line arguments """
    gbff_file: str
    rgi_file: str
    depth_file: str
    sample_name: str
    out_file: str


def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Format rgi file ',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('gbff_file',
                        metavar='bakta gbff output',
                        help='The CDS called by annotator in gbff')

    parser.add_argument('-r',
                        '--rgi',
                        metavar='rgi tsv file',
                        type=str,
                        help='The tsv formatted file produced by rgi',
                        default='',
                        required=True)
    
    parser.add_argument('-d',
                        '--dep',
                        metavar='depth tsv file',
                        type=str,
                        help='The tsv formatted depth file of CDSs produced by the calculate_CDS_depth script',
                        default='',
                        required=True)
    
    parser.add_argument('-n',
                        '--nam',
                        metavar='sample name',
                        type=str,
                        help='The sample name',
                        default='',
                        required=True)

    parser.add_argument('-o',
                        '--out',
                        metavar='output file',
                        type=str,
                        help='name of the output file',
                        default='',
                        required=True)

    args = parser.parse_args()

    if not os.path.isfile(args.gbff_file):
        parser.error("gbff file does not exist")

    if not os.path.isfile(args.rgi):
        parser.error("rgi file does not exist")

    if not os.path.isfile(args.dep):
        parser.error("Depth file does not exist")
    

    return Args(args.gbff_file, args.rgi, args.dep, args.nam, args.out)


args = get_args()

records = SeqIO.parse(args.gbff_file, "genbank") 
table = pandas.read_csv(args.rgi_file, sep="\t", header=0, index_col=None)
sample_name = args.sample_name

# contig	gene	start	end	depth	ratio_assembly	contig_depth	contig_ratio_depth	contig_length
gene_depth_file = pandas.read_csv(args.depth_file, sep="\t", header=0, index_col="gene")

locus2data = {}
for record in records:
    for feature in record.features:
        if feature.type == 'CDS':
            locus2data[feature.qualifiers["locus_tag"][0]] = {"start": int(feature.location.start), 
                                                              "stop": int(feature.location.end), 
                                                              "strand": feature.location.strand, 
                                                              "nucl_seq": str(feature.extract(record.seq)),
                                                              "contig": record.name, 
                                                              "depth": int(gene_depth_file.loc[feature.qualifiers["locus_tag"][0],"depth"]), 
                                                              "contig_depth": int(gene_depth_file.loc[feature.qualifiers["locus_tag"][0],"contig_depth"]),
                                                              "contig_length": int(gene_depth_file.loc[feature.qualifiers["locus_tag"][0],"contig_length"])}


table["Gene_depth"] = None
table["Contig_depth"] = None
table["Contig_length"] = None
table["Sample"] = sample_name
table["reference_db"] = "CARD"

header = {

"ORF_ID": "ORF_ID",
"Contig": "Contig",	
"Start": "Start",
"Stop": "Stop",	
"Orientation": "Strand",	
"Cut_Off": "Cut_Off",
"Pass_Bitscore": "Pass_Bitscore",
"Best_Hit_Bitscore": "Bitscore",
"Best_Hit_ARO": "Best_hit",
"Best_Identities": "Percent_identity",
"ARO": "ARO",
"Model_type": "Model_type",	
"SNPs_in_Best_Hit_ARO": "SNPs",	
"Other_SNPs": "Other_SNPs",
"Drug Class": "Drug_class",
"Resistance Mechanism": "Mechanism",
"AMR Gene Family": "AMR_family",
"Predicted_DNA": "Nucleotide_sequence",
"Predicted_Protein": "Protein_sequence",
"CARD_Protein_Sequence": "CARD_Protein_Sequence",
"Percentage Length of Reference Sequence": "Percent_coverage",
"ID": "ID",
"Model_ID": "Model_ID",	
"Nudged": "Nudged",
"Note": "Note",
	
}

table = table.rename(columns=header)

for n,row in table.iterrows():
    orf_id = row.ORF_ID.split(" ")[0]

    table.loc[n, "ORF_ID"] = orf_id
    table.loc[n, "Start"] = locus2data[orf_id]["start"]
    table.loc[n, "Contig"] = locus2data[orf_id]["contig"]
    table.loc[n, "Stop"] = locus2data[orf_id]["stop"]
    table.loc[n, "Strand"] = locus2data[orf_id]["strand"]
    table.loc[n, "Nucleotide_sequence"] = locus2data[orf_id]["nucl_seq"]
    table.loc[n, "Gene_depth"] = locus2data[orf_id]["depth"]
    table.loc[n, "Contig_depth"] = locus2data[orf_id]["contig_depth"]
    table.loc[n, "Contig_length"] = locus2data[orf_id]["contig_length"]

table["Start"] = table["Start"].astype(int)
table["Stop"] = table["Stop"].astype(int)
table["Strand"] = table["Strand"].astype(int)
table["Gene_depth"] = table["Gene_depth"].astype(int)
table["Contig_depth"] = table["Contig_depth"].astype(int)
table["Contig_length"] = table["Contig_length"].astype(int)


table[["Sample", "reference_db", "ORF_ID", "Contig", "Contig_length", "Contig_depth", "Start", "Stop", "Strand", "Gene_depth","Best_hit", "ARO", "Bitscore", "Pass_Bitscore", "Cut_Off","Percent_identity", "Percent_coverage", "AMR_family", "Model_type", "Mechanism", "SNPs", "Drug_class", "Protein_sequence", "Nucleotide_sequence", "Note"]].to_csv(args.out_file, sep="\t", index=None)

