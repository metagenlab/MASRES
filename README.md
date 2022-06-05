### FlexAMR

FlexAMR is a diagnostic sequencing pipeline that can provide resistance gene predictions from either raw Nanopore reads, Illumina reads or both.

This tool is still under development.

To run this tool, you will need to install the latest version of nextflow (tested on version 22.04.3) and an installation of Singularity (tested on singularity version 3.7.1)

More info on the installation of these tools can be found here:

https://www.nextflow.io/docs/latest/getstarted.html

https://sylabs.io/guides/3.0/user-guide/installation.html

Once those two dependencies have been installed, clone the directory. 

```
git clone -b master https://github.com/metagenlab/MASRES.git

```

After this, you will need to download the database for Homopolish.

```
cd MASRES/
mkdir database
cd database/
wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/bacteria.msh.gz
gunzip bacteria.msh.gz
```
If you've saved the database in a different place, open the config file and edit the path to the database accordingly

Now you have have all you need to start running FlexAMR

FlexAMR has three modes "short" for short read assemblies using Illumina only, "long" for long read assemblies using Oxford Nanopore and
"hybrid" for both Illumina and Oxford Nanopore. Keep in mind that each mode has its own dedicated sample csv format. Some examples of sample csv files 
are provided, all that is required is to replace the names of the samples and the directories of the raw reads.

## Overview of workflow
![Untitled Diagram drawio (4)](https://user-images.githubusercontent.com/70012389/172053345-35fc312a-5ada-4cab-8a62-9aa6ad447d02.png)


## Example of command:

To get help documentation:

nextflow run main.nf --help

To run a hybrid assembly:

nextflow run main.nf --reads_csv samples_hybrid.csv --outdir results_hybrid --threads 10 --mode hybrid

To run short read assembly:

nextflow run main.nf --reads_csv samples_short.csv --outdir results_short --threads 10 --mode short

TO run long read assembly:

nextflow run main.nf --reads_csv samples_long.csv --outdir results_long --threads 10 --mode long
