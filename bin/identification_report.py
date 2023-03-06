#!/opt/conda/bin/python

import argparse
import pandas
import re
import os
import io
from typing import NamedTuple


STYLE = """
        <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap.min.css"/>
        <style type="text/css">
        body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;padding-bottom:10px;background-color:#fff;color:#333;margin:0}body>div .section::before{content:"";display:block;height:80px;margin:-80px 0 0}#identification::before{margin:0}.topic-title{font-size:18pt}body>div>.section{margin-left:255px;margin-bottom:3em}div.section{margin-right:20px}#contents>p{display:none}button,li p.first{display:inline-block}#contents{margin-top:80px;padding-left:0;width:235px;background-color:#f1f1f1;height:100%;position:fixed;overflow:auto}#contents ul{list-style-type:none}#contents ul>li{font-size:14pt}#contents ul>li a:hover{color:#151d26}button,h1.title{color:#fff;background-color:#151d26}#contents ul>li>ul>li{font-size:12pt}h1.title{margin-top:0;position:fixed;z-index:10;padding:20px;width:100%}code,table tr:nth-child(2n),tt{background-color:#f8f8f8}.one-col{min-width:310px;height:500px;margin:0 auto}.two-col-left{height:300px;width:49%;float:left}.two-col-right{height:300px;width:49%;float:right}button{margin:0 5px 0 0;padding:5px 25px;font-size:18px;line-height:1.8;appearance:none;box-shadow:none;border-radius:3px;border:none}button:focus{outline:0}button:hover{background-color:#4183C4}button:active{background-color:#27496d}.legend-rect{width:20px;height:20px;margin-right:8px;margin-left:20px;float:left;-webkit-border-radius:2px;border-radius:2px}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}dl,dl dt,dl dt:first-child,hr,table,table tr{padding:0}table tr td,table tr th{border:1px solid #ccc;text-align:left;padding:6px 13px}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#151d26;border-bottom:1px solid #ccc}h2{font-size:24px;color:#000}h3{font-size:18px}h4{font-size:16px}dl dt,h5,h6{font-size:14px}h6{color:#777}blockquote,dl,li,ol,p,pre,table,ul{margin:15px 0}hr{background:url(http://tinyurl.com/bq5kskr) repeat-x;border:0;color:#ccc;height:4px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}dl dt{font-weight:700;font-style:italic;margin:15px 0 5px}blockquote>:first-child,dl dd>:first-child,dl dt>:first-child,table tr td :first-child,table tr th :first-child{margin-top:0}blockquote>:last-child,dl dd>:last-child,dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}table{border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0}table tr th{font-weight:700;margin:0}table tr td{margin:0}table tr td :last-child,table tr th :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame>span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center>span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right>span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right>span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;background:0 0}.highlight pre,pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}h1{line-height:1.6}.simple{padding-left:20px}.docutils.container{width:100%}
        .pull-left{
        .dataTables_filter {
        float: left !important;
        }
        </style>
    """

SCRIPT = """
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
        <script src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/1.10.19/js/dataTables.bootstrap.min.js"></script>
        <script src="https://unpkg.com/cytoscape/dist/cytoscape.min.js"></script>
        <script src="https://unpkg.com/webcola/WebCola/cola.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/cytoscape-cola@2.2.4/cytoscape-cola.min.js"></script>

        <script>

        $(document).ready(function() {
            $('#gtdbtk_table').DataTable( {
                dom: 'Bfrtip',
                "pageLength": 10,
                "order": [[1, 'desc']],
                "searching": true,
                "bLengthChange": false,
                "paging":   true,
                "info": false
            } );
        } );

        $(document).ready(function() {
            $('#centri_table').DataTable( {
                dom: 'Bfrtip',
                "pageLength": 10,
                "order": [[0, 'desc']],
                "searching": true,
                "bLengthChange": false,
                "paging":   true,
                "info": false
            } );
        } );
        
        $(document).ready(function() {
            $('#checkm_table').DataTable( {
                dom: 'Bfrtip',
                "pageLength": 10,
                "order": [[0, 'desc']],
                "searching": true,
                "bLengthChange": false,
                "paging":   true,
                "info": false
            } );
        } );
                
        $(document).ready(function() {
            $('#mash_table').DataTable( {
                dom: 'Bfrtip',
                "order": [[1, 'desc']],
                "pageLength": 10,
                "searching": true,
                "bLengthChange": false,
                "paging":   true,
                "info": false,
                'rowCallback': function(row, data, index){
                    $(row).find('td:eq(2)').css('background-color', 'rgba(255, 0, 0, 0.2)');
                    $(row).find('td:eq(3)').css('background-color', 'rgba(255, 0, 0, 0.2)');
                    $(row).find('td:eq(4)').css('background-color', 'rgba(255, 0, 0, 0.2)');
                    $(row).find('td:eq(5)').css('background-color', 'rgba(0,128,0, 0.2)');
                    $(row).find('td:eq(6)').css('background-color', 'rgba(0,128,0, 0.2)');
                    $(row).find('td:eq(7)').css('background-color', 'rgba(0,128,0, 0.2)');
                    $(row).find('td:eq(8)').css('background-color', 'rgba(128,128,128, 0.2)');
                    $(row).find('td:eq(9)').css('background-color', 'rgba(128,128,128, 0.2)');
                    $(row).find('td:eq(10)').css('background-color', 'rgba(128,128,128, 0.2)');
                },
            } );
        } );
        </script>

    """


class Report():
    
    def __init__(self, report_type="Identification"):        
        
        
        self.summary_table = []
        
        
        self.report_str = f"""

.. raw:: html

    {SCRIPT}

    {STYLE}


Diag Pipeline - {report_type} report
=============================================================

.. contents::
    :backlinks: none
    :depth: 2

        """
    
    
    def add_section(self, section_name):
        underl = '-' * len(section_name)
        
        self.report_str += f'''
{section_name}
{underl}
        '''
        
        
    def add_subsection(self, 
                       subsection_name):
        underl = '*' * len(subsection_name)
        
        self.report_str += f'''
{subsection_name}
{underl}
        '''
        
    def add_table(self, table_html):
        
        
        self.report_str += f'''
.. raw:: html

    {table_html}
        '''
        
    def checkm_table(self, checkm_tables):
        
        df_list = [pandas.read_csv(checkm_table, delimiter="\t", header=0) for checkm_table in checkm_tables]
        
        df = pandas.concat(df_list)
        df = df.rename(columns={"Bin Id":"Sample"})
        df = df[["Sample","Marker lineage","# markers", "Completeness","Contamination","Strain heterogeneity","0","1","2","3","4","5+"]]
        pandas.set_option('display.max_colwidth', None)


        # rename columns for standardized report and upload to LIMS
        df["2+"] = df[["2","3", "4", "5+"]].sum(axis=1)
        cols = {
            "Sample": "sample",
            "Marker lineage":"checkm_lineage",
            "# markers":"checkm_markers", 
            "Completeness":"checkm_completeness",
            "Contamination":"checkm_contamination",
            "Strain heterogeneity":"checkm_heterogeneity",
            "0":"checkm_markers_missing",
            "1":"checkm_markers_single",
            "2+":"checkm_markers_multiple",
        }

        df_resh = df[cols.keys()].rename(columns=cols).set_index(["sample"]).stack()
        self.summary_table += df_resh.reset_index().values.tolist()


        df = df.drop("2+", axis=1)
        
        df_str = df.to_html(
            index=False,
            bold_rows=False,
            classes=["dataTable"],
            table_id="checkm_table",
            escape=False,
            border=0)

        return df_str.replace("\n", "\n" + 10 * " ")
    

    def get_mash_table(self, 
                       file_list, 
                       keep=1):
        '''
        sample
        hit 1: score (shared): name
        hit 2:
        hit 3:

        '''
        df_list = []
        for one_file in file_list:
            df = pandas.read_csv(one_file)[["query_name", "intersect_bp","f_orig_query","f_match","name","f_match_orig"]]
            df = df.rename(columns={"query_name":"Sample"})
            df_list.append(df)
        df_merged = pandas.concat(df_list)
        df_merged["Sample"] = df_merged["Sample"].astype(str)
        df_merged["Sample"] = [i.split(".fna")[0] for i in df_merged["Sample"]]
        df_merged["f_orig_query"] = [round(i * 100, 2) for i in df_merged["f_orig_query"]]
        df_merged["f_match"] = [round(i * 100, 2) for i in df_merged["f_match"]]
        df_merged["acc"] = [f'<a href="https://www.ncbi.nlm.nih.gov/assembly/{i.split(" ")[0]}">{i.split(" ")[0]}</a>' for i in df_merged["name"]]
        df_merged["description"] = [' '.join(i.split(" ")[1:]) for i in df_merged["name"]]
        df_merged = df_merged.drop(["name", "f_match_orig"], axis=1)
        df_merged_filt = df_merged.groupby(["Sample"]).head(1).reset_index()

        pandas.set_option('display.max_colwidth', None)

        df_str = df_merged_filt.reset_index(drop=True).to_html(
            index=False,
            bold_rows=False,
            classes=["dataTable"],
            table_id="mash_table",
            escape=False,
            border=0)

        return df_str.replace("\n", "\n" + 10 * " ")


    def get_skani_table(self, 
                        skani_tsv, 
                        keep=1):
        '''
        sample
        hit 1: score (shared): name
        hit 2:
        hit 3:

        '''
        df = pandas.read_csv(skani_tsv, header=0, sep="\t")
        df = df[["Query_file", "Ref_file", "ANI", "Align_fraction_ref", "Align_fraction_query", "Ref_name"]]
        df = df.rename(columns={"Query_file":"Sample", "Align_fraction_ref": "fraction ref", "Align_fraction_query": "fraction query"})
        # Ref_file	Query_file	ANI	Align_fraction_ref	Align_fraction_query	Ref_name	Query_name
        
        df["Sample"] = [i.split(".fna")[0] for i in df["Sample"]]
        df["Ref_file"] = [i.split("/")[-1].split("_genomic.fna.gz")[0] for i in df["Ref_file"]]
        #df["Ref_name"] = [i.split(".fna")[0] for i in df["Ref_name"]]
        #df["f_orig_query"] = [round(i * 100, 2) for i in df["Align_fraction_query"]]
        #df["f_match"] = [round(i * 100, 2) for i in df["Align_fraction_ref"]]
        
        df_filt = df.groupby(["Sample"]).head(keep).reset_index()
        df_filt["acc"] = [f'<a href="https://www.ncbi.nlm.nih.gov/assembly/{i}">{i}</a>' for i in df_filt["Ref_file"]]
        #df_merged_filt["acc"] = [f'<a href="https://www.ncbi.nlm.nih.gov/nucleotide/{i.split(" ")[0]}">{i.split(" ")[0]}</a>' for i in df_merged_filt["Ref_name"]]
        df_filt["description"] = [' '.join(i.split(" ")[1:]) for i in df_filt["Ref_name"]]
        
        df_filt = df_filt.drop("Ref_name", axis=1)


        # extract columns for tsv report
        ###############
        cols = {
            "Sample": "sample",
            "fraction ref":"skani_fraction_ref",
            "fraction query":"skani_fraction_query", 
            "ANI":"skani_ANI",
            "description":"skani_description",
        }

        df_resh = df_filt[cols.keys()].rename(columns=cols).set_index(["sample"]).stack()

        self.summary_table += df_resh.reset_index().values.tolist()


        pandas.set_option('display.max_colwidth', None)

        df_str = df_filt.reset_index(drop=True).to_html(
            index=False,
            bold_rows=False,
            classes=["dataTable"],
            table_id="skani_table",
            escape=False,
            border=0)

        return df_str.replace("\n", "\n" + 10 * " ")


    def get_multiqc_table(self, 
                          assembly_multiqc=False,
                          mapping_multiqc=False):

        mq_table = []
        if assembly_multiqc:
            mq_table.append(["MultiQC genome assemblie(s)", '<a href="%s">MultiQC</a>' % '/'.join(assembly_multiqc.split('/')[1:])])

        if mapping_multiqc:
            for n, multiqc in enumerate(mapping_multiqc):
                multiqc_link = '<a href="%s">MultiQC</a>' % '/'.join(multiqc.split('/')[1:])
                mq_table.append(["%s" % re.sub("_", " ", multiqc.split("/")[1]), multiqc_link])
        header = ["Name", "Link"]

        df = pandas.DataFrame(mq_table, columns=header)

        # cell content is truncated if colwidth not set to -1
        pandas.set_option('display.max_colwidth', None)

        df_str = df.to_html(
            index=False,
            bold_rows=False,
            classes=["dataTable"],
            table_id="multiqc_table",
            escape=False,
            border=0)

        return df_str.replace("\n", "\n" + 10 * " ")
    

    def get_rrna_summary_table(self, 
                            raw_table):
        
        df = pandas.read_csv(raw_table, delimiter="\t", header=0)
        
        # sample	query	hit	alignment_length	alignment_length	percent_identity	evalue	bitscore
        row_list = []
        sample2count = {}
        for n, row in df.iterrows():
            sample = str(row["sample"])
            contig = row["query"]
            BBH_taxnonomy = row["hit"].split(",")[1:]
            identity = row["percent_identity"]
            if sample not in sample2count:
                sample2count[sample] = 0
            sample2count[sample] += 1
            row_list.append([sample,
                            sample2count[sample], 
                            contig, 
                            BBH_taxnonomy, 
                            identity])
        
        header = ["sample", "N.", "expected", "contig", "BBH_taxnonomy", "identity"]
        df = pandas.DataFrame(row_list, columns=header)
        pandas.set_option('display.max_colwidth', -1)

        df_str = df.to_html(
            index=False,
            bold_rows=False,
            classes=["dataTable"],
            table_id="16S_table",
            escape=False,
            border=0)

        return df_str.replace("\n", "\n" + 10 * " ")
    
    def parse_gtdbtk(self, gtdbtk_summary):
        
        df = pandas.read_csv(gtdbtk_summary, header=0, sep="\t")
        
        df = df[["user_genome", "classification","closest_placement_taxonomy", "closest_placement_ani", "closest_placement_af", "closest_placement_reference", "classification_method"]]
        df = df.rename(columns={"user_genome":"Sample", "closest_placement_ani":"ani", "closest_placement_af":"alignment fraction", "closest_placement_taxonomy":"Closest", "closest_placement_reference": "Closest acc."})
        df["Closest"] = [i.split(";")[-1] for i in df["Closest"]]
        df["acc."] = [f'<a href="https://www.ncbi.nlm.nih.gov/assembly/{i}">{i}</a>' for i in df["Closest acc."]]
        df = df.rename(columns={"user_genome":"Sample"})
        df["classification"] = [', '.join(i.split(";")) for i in df["classification"]]
        
        # extract columns for tsv report
        ###############
        cols = {
            "Sample": "sample",
            "ani":"gtdbtk_closest_ani",
            "alignment fraction":"gtdbtk_closest_fraction_ref", 
            "Closest":"gtdbtk_closest_taxonomy",
            "Closest acc.":"gtdbtk_closest_accession",
            "classification_method": "gtdbtk_classification_method",
            "classification": "gtdbtk_classification"
        }

        df_resh = df[cols.keys()].rename(columns=cols).set_index(["sample"]).stack()
        self.summary_table += df_resh.reset_index().values.tolist()
        ###############
        
        # generate html page
        df_str = df.to_html(
            index=False,
            bold_rows=False,
            classes=["dataTable"],
            table_id="gtdbtk_table",
            escape=False,
            border=0)

        return df_str.replace("\n", "\n" + 10 * " ")   
    
    def get_centrifuge_table(self, centrifuge_links):
        row_list = []

        if not os.path.exists("centrifuge"):
            os.makedirs("centrifuge")
        # N2405.reportk.txt
        # # ["percent", "n_reads", "n_reads_direct", "rank", "taxid", "scientific_name"]
        for sample in centrifuge_links:
            sample_name = sample.split('.reportk.txt')[0]
            
            table = pandas.read_csv(sample, delimiter="\t", names=["percent_root", "number_rooted", "number_assigned","rank", "taxid", "scientific_name"])
           
            total_alligned_reads = sum(table["number_assigned"])
            table["percent_assigned"] = round((table["number_assigned"] / total_alligned_reads) * 100, 2)
            
            table.sort_values(["percent_assigned"], ascending=False).to_csv(f"centrifuge/{sample}_centrifuge.tsv", sep="\t")
            table_filter = table.sort_values(["percent_assigned"], ascending=False).query('rank=="G" or rank=="S"')
            #unclassified = table.iloc[0]
            hit_1 = table_filter.iloc[0]
            hit_2 = table_filter.iloc[1]
            hit_3 = table_filter.iloc[2]

            row = [sample_name,
                   "n/a",
                    "%s (%s)" % (hit_1["scientific_name"], hit_1["rank"]),
                    hit_1["percent_assigned"],
                    "%s (%s)" % (hit_2["scientific_name"], hit_2["rank"]),
                    hit_2["percent_assigned"],
                    "%s (%s)" % (hit_3["scientific_name"], hit_3["rank"]),
                    hit_3["percent_assigned"],
                    f'<a href="centrifuge/{sample}_centrifuge.tsv">detail</a>']

            row_list.append(row)
        header = ["Sample", "expected", "taxon 1", "% reads", "taxon 2", "% reads", "taxon 3", "% reads", "detail"]
        df = pandas.DataFrame(row_list, columns=header).reset_index(drop=True)
        pandas.set_option('display.max_colwidth', None)

        df_str = df.to_html(
            index=False,
            bold_rows=False,
            classes=["dataTable"],
            table_id="centri_table",
            escape=False,
            border=0)

        return df_str.replace("\n", "\n" + 10 * " ")
    

    def publish(self, out_prefix):
        from docutils.core import publish_file
        with open(f"{out_prefix}.rst", "w") as fr:
            fr.write(self.report_str)
        with open(f"{out_prefix}.html", "w") as fh:
            publish_file(
                source=io.StringIO(self.report_str),
                destination=fh,
                writer_name="html",
                settings_overrides={"stylesheet_path": ""},
            )

    def write_tsv(self, out_prefix):
        df = pandas.DataFrame(self.summary_table, columns=["sample", "term", "value"])
        df.to_csv(f"{out_prefix}.tsv", sep="\t", index=False)
    
    
if __name__ == "__main__":
    
    from pathlib import Path
    import argparse
    
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate report",
        epilog="",
    )
    parser.add_argument(
        "-m", '--mash',
        nargs="+",
        help="SourMash files",
    )
    parser.add_argument(
        "-c", '--centrifuge',
        nargs="+",
        help="Centrifuge files.",
    )
    parser.add_argument(
        "-g", '--gtdbtk_summary',
        type=Path,
        help="GTDBtk summary file",
    )
    parser.add_argument(
        "-k", '--checkm',
        nargs="+",
        help="CheckM results",
    )
    parser.add_argument(
        "-s", '--skani',
        help="skani tsv",
    )
    parser.add_argument(
        "-o", '--out',
        type=str,
        help="Output prefix",
    )



    args = parser.parse_args()

    report = Report()
    
    report.add_section("identification")

    # gtdbtk_summary
    report.add_subsection("GTDBtk Classification (assembly)")
    report.add_table(report.parse_gtdbtk(args.gtdbtk_summary))

    # mash
    report.add_subsection("SourMash vs GTDB (assembly)")
    report.add_table(report.get_mash_table(args.mash))
 
     # skani
    report.add_subsection("Skani vs GTDB (assembly)")
    report.add_table(report.get_skani_table(args.skani))
    
    # centrifuge
    report.add_subsection("Centrifuge classification (reads; Genus and Species level)")
    report.add_table(report.get_centrifuge_table(args.centrifuge))
    
    # checkm 
    report.add_subsection("CheckM markers QC")
    report.add_table(report.checkm_table(args.checkm))
    
    # print
    report.publish(args.out)
    report.write_tsv(args.out)