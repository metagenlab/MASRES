#!/bin/bash

### get bakta db
mkdir bakta_db
wget https://zenodo.org/record/7025248/files/db.tar.gz?download=1 -P bakta_db

### get homopolish db
mkdir homopolish_db
wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/bacteria.msh.gz -P homopolish_db

### get platon db
mkdir platon_db
wget https://zenodo.org/record/4066768/files/db.tar.gz -P platon_db

### get sourmash db
mkdir sourmash_db
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic.k21.lca.json.gz -P sourmash_db

