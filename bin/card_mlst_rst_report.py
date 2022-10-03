#!/opt/conda/envs/python-r/bin/python
import argparse
import pandas
import re
import os
import io
from typing import NamedTuple

class Args(NamedTuple):
    """ Command-line arguments """
    mlst_file: str
    rgi_file: str
    mash_file: str
    plasmid_file: str
    out_file: str


def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Format merged report ',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('mlst_file',
                        metavar='mlst file',
                        help='The mlst tsv file produced by MLST')

    parser.add_argument('-r',
                        '--rgi',
                        metavar='rgi tsv files',
                        type=str,
                        help='The tsv formatted file produced by rgi',
                        default='',
			nargs='+',
                        required=True)
    
    parser.add_argument('-m',
                        '--mas',
                        metavar='mash files',
                        type=str,
                        help='The tsv formatted file produced by MASH SCREEN',
                        default='',
			nargs='+',
                        required=True)
    
    parser.add_argument('-p',
                        '--pla',
                        metavar='plasmid files',
                        type=str,
                        help='The predicted plasmids file produced by PLATON',
                        default='',
			nargs='+',
                        required=True)

    parser.add_argument('-o',
                        '--out',
                        metavar='output file',
                        type=str,
                        help='name of the output file',
                        default='combined_report.rst',
                        )

    args = parser.parse_args()


    return Args(args.mlst_file, args.rgi, args.mas, args.pla, args.out)





def parse_mlst(mlst_tsv_file):
    sample2mlst = {}
    with open(mlst_tsv_file, "r") as f:
        for row in f:
            data = row.rstrip().split("\t")
            # samples/N1097/annotation/N1097.fna
            sample = data[0]
            mlst_scheme = data[1]
            mlst = data[2]
            sample2mlst[sample] = [mlst_scheme, mlst]
    return sample2mlst


def parse_rgi(rgi_file_list, plasmid_contig_list,
              query_cov_cutoff=30):

    BETALACTAMS = ["monobactam", "carbapenem", "penam", "cephem", "penem", "cephamycin", "cephalosporin", "Beta-Lactam"]
    # carbapenem; cephalosporin; cephamycin; penam
    # create dictionary with sample name as key and list of contigs classified as being on plasmids as values
    sample2plasmid = {}
    for plasmid_file in plasmid_contig_list:
        samplebase = os.path.basename(plasmid_file)
        sample = os.path.splitext(samplebase)[0].split("-")[1]
        t = pandas.read_csv(plasmid_file, sep="\t", header=0)
        contig_plasmids = list(t["ID"])
        sample2plasmid[sample] = contig_plasmids

    sample2rgi = {}
    for rgi_file in rgi_file_list:
        print("rgi_file", rgi_file)
        samplebase = os.path.basename(rgi_file)
        sample = os.path.splitext(samplebase)[0]
        plasmids = sample2plasmid[sample]
        sample2rgi[sample] = {}
        sample2rgi[sample]["transporters"] = []
        sample2rgi[sample]["SNP"] = []
        sample2rgi[sample]["drug_resistance"] = {}
        t = pandas.read_csv(rgi_file, sep="\t", header=0)
        for n, row in t.iterrows():
            # print("-------", row)
            if row["Contig"] in plasmids:
                gene = f"{row['Best_hit']} (p)"
            else:
                gene = f"{row['Best_hit']} (c)"

            model_type = row["Model_type"]
            mechanism = row["Mechanism"]
            if not isinstance(mechanism, str):
                mechanism = 'n/a'
            coverage = float(row["Percent_coverage"])
            if float(coverage) < query_cov_cutoff:
                print("Skipping low cov entry: %s (%s %%)" % (gene, coverage))
                continue
            identity = float(row["Percent_identity"])
            if "efflux" in mechanism:
                sample2rgi[sample]["transporters"].append([gene, coverage, identity])
                continue

            # protein variant model
            # protein homolog model
            # rRNA gene variant model
            if model_type in ["protein variant model", "rRNA gene variant model"]:
                SNPs_in_Best_Hit_ARO = row["SNPs"]
                sample2rgi[sample]["SNP"].append([gene, SNPs_in_Best_Hit_ARO, coverage, identity])
                continue
            elif model_type == "protein homolog model":
                try:
                    drug_class_list = row["Drug_class"].split("; ")
                except:
                    drug_class_list = ["Unspecified"]
                drug_class_list = list(set([i if i not in BETALACTAMS else "Beta-Lactam" for i in drug_class_list]))

                for drug in drug_class_list:
                    drug = re.sub(" antibiotic", "", drug)
                    if drug in BETALACTAMS:
                        drug = "Beta-Lactam"
                    if drug not in sample2rgi[sample]["drug_resistance"]:
                        sample2rgi[sample]["drug_resistance"][drug] = {}
                    if gene not in sample2rgi[sample]["drug_resistance"][drug]:
                        sample2rgi[sample]["drug_resistance"][drug][gene] = [1, coverage, identity]
                    else:
                        # multiple copies of the same gene, increment count and keep lowest coverage and identity
                        sample2rgi[sample]["drug_resistance"][drug][gene][0] += 1
                        # keep lowest coverage
                        if coverage < sample2rgi[sample]["drug_resistance"][drug][gene][1]:
                            sample2rgi[sample]["drug_resistance"][drug][gene][1] = coverage
                        # keep lowest identity
                        if coverage < sample2rgi[sample]["drug_resistance"][drug][gene][2]:
                            sample2rgi[sample]["drug_resistance"][drug][gene][2] = identity
            else:
                #
                raise IOError("Unknown model type:", model_type)
    return sample2rgi


def parse_mash(mash_file_list):
    sample2species = {}
    for mash_file in mash_file_list:
        samplebase = os.path.basename(mash_file)
        sample = os.path.splitext(samplebase)[0]
        # 0.999905	998/1000	0	Klebsiella pneumoniae strain 196 map unlocalized plasmid unnamed1 Plasmid_1_Contig_2, whole genome shotgun sequence
        with open(mash_file, "r") as f:
            rows = [i.rstrip().split(",") for i in f]
            best_hit_score = rows[1][1]
            best_hit_description = ' '.join(rows[1][9].split(" ")[1:3])
            sample2species[sample] = [best_hit_description, best_hit_score]
    return sample2species
    



def generate_report_rst(sample2mlst, 
                        sample2rgi, 
                        sample2species,
                        output_name):

    print("sample2species within", sample2species)
    ######################
    # MLST table
    ######################
    SAMPLES_LIST = list(sample2species.keys())
    MASH_DATA = []
    for sample in SAMPLES_LIST:
        species = sample2species[sample][0]
        score = sample2species[sample][1]
        MASH_DATA.append(f"{species} ({score})")
    MLST_list = [f"{scheme},{mlst}" for scheme, mlst in sample2mlst.values()]

    table_1_rows = [','.join(i) for i in zip(SAMPLES_LIST, MLST_list, MASH_DATA)]
    table_1 = '\n    '.join(table_1_rows)


    ######################
    # Transporters table
    ######################
    transporters_table_rows = []
    for sample in SAMPLES_LIST:
        if sample in sample2rgi:
            if len(sample2rgi[sample]["transporters"]) != 0:
                transporter_str_list = []
                for tranporter, coverage, identity in sample2rgi[sample]["transporters"]:
                    if identity < 90:
                        transporter_str = f"{tranporter}[1]"
                    else:
                        transporter_str = f"{tranporter}"
                    if coverage < 80:
                        transporter_str += f"[2]"
                    transporter_str_list.append(transporter_str)
                transporters = '; '.join(transporter_str_list)
            else:
                transporters = "No known transporters found"
            transporters_table_rows.append(f"* - {sample}\n      - {transporters}")
        transporters_table = '\n    '.join(transporters_table_rows)


    ######################
    # SNP table
    ######################

    SNP_table_rows = []
    for sample in SAMPLES_LIST:
        if sample in sample2rgi:
            if len(sample2rgi[sample]["SNP"]) != 0:
                SNPs_str_list = []
                for gene, snp, coverage, identity in sample2rgi[sample]["SNP"]:
                    if identity < 90:
                        SNPs_str = f"{gene} ({snp})[1]"
                    else:
                        SNPs_str = f"{gene} ({snp})"
                    if coverage < 80:
                        SNPs_str += f"[2]"
                    SNPs_str_list.append(SNPs_str)
                SNPs = '\n        '.join(SNPs_str_list)
            else:
                SNPs = "No resistance associated SNPs"
            SNP_table_rows.append(f"* - {sample}\n      - {SNPs}")
        SNP_table = '\n    '.join(SNP_table_rows)

    ######################
    # Resistance table
    ######################

    resistance_table_rows = []
    for sample in SAMPLES_LIST:
        if sample in sample2rgi:
            sample_str = ''
            if len(sample2rgi[sample]["drug_resistance"]) != 0: 
                drug_list = list(sample2rgi[sample]["drug_resistance"].keys())
                drug_list.sort(key=lambda v: v.upper())
                for drug in drug_list:
                    gene_str_list = []
                    drug_format = f"{drug[0].upper()}{drug[1:len(drug)]}"
                    for gene in sample2rgi[sample]["drug_resistance"][drug]:
                        gene_data = sample2rgi[sample]["drug_resistance"][drug][gene]
                        gene_count, coverage, identity = gene_data
                        # indicate if multiple copies of the same gene
                        if gene_count > 1:
                            gene_label = f'{gene} ({gene_count}x)'
                        else:
                            gene_label = gene
                        if identity < 90:
                            gene_str = f"{gene_label}[1]"
                        else:
                            gene_str = f"{gene_label}"
                        if coverage < 80:
                            gene_str += f"[2]"
                        gene_str_list.append(gene_str)
                    sample_str += f' | **{drug_format}:** ' + ", ".join(gene_str_list) + '\n        '
            else:
                gene_str = "No resistance genes found"
            resistance_table_rows.append(f"* - {sample}\n      - {sample_str}")
        resistance_table = '\n    '.join(resistance_table_rows)
    report_str = f"""

==================
Resistance report
==================

Strain identification 
---------------------

.. csv-table::
    :header: "Sample", "Scheme", "MLST", "Mash best hit"
    :widths: 7, 10, 7, 30

    {table_1}


Antibiotic resistance genes
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 10, 50
   :header-rows: 1

    * - **Sample**
      - **Antibiotic resistance genes**
    {resistance_table}
    
| [1] identity < 90%
| [2] partial gene (<80% of the length of the reference)
| [c] gene present on the chromosome
| [p] gene present on a plasmid

Single nucleotide polymorphisms (SNP)
--------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 10, 50
   :header-rows: 1

    * - **Sample**
      - **SNP(s)**
    {SNP_table}

| [1] identity < 90%
| [2] partial gene (<80% of the length of the reference)
 
Antibiotic efflux systems (& regulators)
-----------------------------------------

.. list-table::
    :header-rows: 1
    :widths: 10, 50
    
    * - **Sample**
      - **Transporters**
    {transporters_table}

| [1] identity < 90%
| [2] partial gene (<80% of the length of the reference)
| [c] gene present on the chromosome
| [p] gene present on a plasmid

    """



    with open(output_name, "w") as fh:
        fh.write(report_str)

args = get_args()

rgi_file_list = args.rgi_file
plasmid_file_list = args.plasmid_file
# input either list of file or a single file
if not isinstance(rgi_file_list, list):
    rgi_file_list = [rgi_file_list]
mlst_file = args.mlst_file
mash_file_list = args.mash_file
# input either list of file or a single file
if not isinstance(mash_file_list, list):
    mash_file_list = [mash_file_list]

output_file = args.out_file
   
sample2mlst = parse_mlst(mlst_file)
sample2rgi = parse_rgi(rgi_file_list, plasmid_file_list)

sample2species = parse_mash(mash_file_list)

sample2mlst = {key:value for key,value in sample2mlst.items() if key in sample2rgi}

generate_report_rst(sample2mlst, 
                    sample2rgi, 
                    sample2species,
                    output_file)
