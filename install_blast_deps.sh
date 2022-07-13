#!/bin/bash
CWD="$(pwd)"
DB_DIR=$CWD/db
SRC_DIR=$CWD/src

echo "Downloading BLAST+ executables and swissprot database..."

#======================================================================#
echo "Downloading BLAST+..."
cd $SRC_DIR
if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform
    BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-macosx.tar.gz"
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Do something under GNU/Linux platform
    BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz"

wget $BLAST_URL
tar -xzf ncbi-blast-2.13.0+-x64-*.tar.gz
rm ncbi-blast-2.13.0+-x64-*.tar.gz
#======================================================================#
echo "Downloading SwissProt DB..."
cd $DB_DIR
wget \
  https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz \
  -O swissprot.tar.gz

tar -zxf swissprot.tar.gz
rm swissprot.tar.gz
#======================================================================#
cd $CWD

