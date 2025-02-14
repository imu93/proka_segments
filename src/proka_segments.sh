#!/bin/bash

SRCPATH=~/storage/Data/cEle_wagos_eco/src/proka_segments/src
GENOME=$1
ID=$(echo $GENOME | sed 's/\..*//')
ANNFIL=${ID}.annotations.gff
GENE=${ID}.disjon_gene.gff3
NCRA=${ID}.disjoin_ncRNA.gff3
TEAN=${ID}.disjoin_repeats.gff3
OUTAN=${ID}.segmented_annotation.gff3

if [[ -s "$GENOME" ]]; then
    echo "Processing '$GENOME'."
else
    echo "The genome file does not exist or is empty."
fi


if [[ -s "$ANNFIL" ]]; then
    echo "Processing '$ANNFIL'."
else
    echo "The genome file does not exist or is empty."
fi

# Disjoin the gene annotation
Rscript $SRCPATH/01.disjoin_protein_coding_genes.R -f $GENOME -g $ANNFIL -o $GENE

# Disjoin the ncRNA annotation
Rscript $SRCPATH/02.reduce_ncrna_genes.R -f $GENOME -g $ANNFIL -o $NCRA

# Disjoin repeat annotation
Rscript $SRCPATH/03.reduce_repeats.R -f $GENOME -g $ANNFIL -o $TEAN

# Now run remove overlaps between all features
Rscript $SRCPATH/04.segemeted_prokka.R -f $GENOME -p $GENE -r $TEAN -n $NCRA -o $OUTAN

rm -rf ${GENOME}.fai

