#!/bin/bash
CWD="$(pwd)"
DB_DIR=${CWD}/db
SRC_DIR=${CWD}/src
OS=$(uname -s)

echo "Downloading BLAST+ executables and swissprot database..."

#======================================================================#
echo "Downloading BLAST+..."
cd "${SRC_DIR}" || exit

if [ "${OS}" == "Darwin" ]; then
	# Do something under Mac OS X platform
	BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-macosx.tar.gz"
elif [ "${OS}" == "Linux" ]; then
	# Do something under GNU/Linux platform
	BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz"
fi

wget "$BLAST_URL" >/dev/null 2>&1

tar -xzf ncbi-blast-2.13.0+-x64-*.tar.gz
rm ncbi-blast-2.13.0+-x64-*.tar.gz
#======================================================================#
echo "Downloading SwissProt DB..."
cd "$DB_DIR" || exit

wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz >/dev/null 2>&1

tar -zxf swissprot.tar.gz
rm swissprot.tar.gz
#======================================================================#
cd "$CWD" || exit
