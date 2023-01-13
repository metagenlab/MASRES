
### MASRES

MASRES is a diagnostic sequencing pipeline that can provide resistance gene predictions from either raw Nanopore reads, Illumina reads or both.

This tool is still under development.

To run this tool, you will need to install the latest version of nextflow (tested on version 22.04.3) and an installation of Singularity (tested on singularity version 3.7.1)

More info on the installation of these tools can be found here:

https://www.nextflow.io/docs/latest/getstarted.html

https://sylabs.io/guides/3.0/user-guide/installation.html

Once those two dependencies have been installed, clone the directory. 

```
git clone -b master https://github.com/metagenlab/MASRES.git

```

After this, you will need to download the databases. A bash script is provided that will download some of these databases:

```
cd MASRES/
bash ./db_download.sh
```
The centrifuge database needs to be downloaded manually using the info on https://ccb.jhu.edu/software/centrifuge/manual.shtml#nt-database

For the BLDB database, since it is currently not available online, the fasta file in the config can be replaced with any custom fasta file with protein sequences of resistance genes.

The config file needs to be edited with the paths to the downloaded databases

Now you have have all you need to start running MASRES

MASRES has three modes "short" for short read assemblies using Illumina only, "long" for long read assemblies using Oxford Nanopore and
"hybrid" for both Illumina and Oxford Nanopore. Keep in mind that each mode expects specific inputs in the sample csv. An example file is provided, all that is required is to replace the names of the samples and the directories of the raw reads and add the corresponding flag when launching the pipeline.

## Overview of workflow
![Untitled Diagram drawio (4)](https://user-images.githubusercontent.com/70012389/212329086-aad5f0bd-b62c-4b0e-8122-1e8115d708c0.png)


## Example of command:

To get help documentation:

```
nextflow run main.nf --help
```

To run a hybrid assembly:

```
nextflow run main.nf --reads_csv samples_hybrid.csv --outdir results_hybrid --mode hybrid
```

To run short read assembly:

```
nextflow run main.nf --reads_csv samples_short.csv --outdir results_short --mode short
```
TO run long read assembly:

```
nextflow run main.nf --reads_csv samples_long.csv --outdir results_long --mode long
```
