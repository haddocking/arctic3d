#!/bin/bash
CWD="$(pwd)"
echo "Current working directory: ${CWD}"
DB_DIR=${CWD}/db
SRC_DIR=${CWD}/src
OS=$(uname -s)
echo "Operating System: ${OS}"
echo "Downloading BLAST+ executables and swissprot database..."

#======================================================================#
echo "Downloading BLAST+..."
cd "${SRC_DIR}" || exit

# if the OS is Darwin or Mac OS X, download the Mac OS X version of BLAST+
# if the OS is Linux, download the Linux version of BLAST+
if [ "${OS}" == "Darwin" ] || [ "${OS}" == "Mac OS X" ]; then
	# Do something under Mac OS X platform
	echo "Downloading BLAST+ for Mac OS X..."
	BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-macosx.tar.gz"
elif [ "${OS}" == "Linux" ]; then
	# Do something under GNU/Linux platform
	echo "Downloading BLAST+ for Linux..."
	BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz"
fi

wget "$BLAST_URL" >/dev/null 2>&1

tar -xzf ncbi-blast-2.15.0+-x64-*.tar.gz
rm ncbi-blast-2.15.0+-x64-*.tar.gz
#======================================================================#
echo "Downloading SwissProt DB..."
cd "$DB_DIR" || exit

wget https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz >/dev/null 2>&1

tar -zxf swissprot.tar.gz
rm swissprot.tar.gz
#======================================================================#
cd "$CWD" || exit
