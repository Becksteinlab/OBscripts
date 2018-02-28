#!/bin/bash
# written 2012 Oliver Beckstein <oliver.beckstein@asu.edu>
# released into the public domain

PROGRAM=$(basename $0)

PDBURL="https://files.rcsb.org/download/"
usage="usage: ${PROGNAME} [options] PDBID [PDBID ...]

Download files PDBID from ${PDBURL}. Note that file names will
have the pdbid in lower case.

Options:

-h             show this help
-d DIRECTORY   place pdbs in this directory
"""

CACHEDIR="."

function die () {
    # die errormessage [error_number]
    local errmsg="$1" errcode="${2:-1}"
    echo "ERROR: ${errmsg}"
    exit ${errcode}
}

function warn () {
    # warn message
    local msg="$1"
    echo "WARNING: ${msg}"
    return
}


#------------------------------------------------------------
# main script
#------------------------------------------------------------

#------------------------------------------------------------
# process options and arguments
#
while getopts hd: OPT; do
    case "${OPT}" in
	h)  echo "${usage}";
	    exit 0
	    ;;
	d)  CACHEDIR="${OPTARG}"
	    ;;
	?)  die "unknown option or missing arument; see -h for usage" 2
	    ;;
    esac
done

shift $((OPTIND - 1))
PDBIDS="$*"
#
#------------------------------------------------------------


#------------------------------------------------------------
# check input
if [ -z "${PDBIDS}" ]; then
    die "No PDBID provided. See '${PROGRAM} -h' for help."
fi

if [ ! -e "${CACHEDIR}" ]; then
    warn "download directory ${CACHEDIR} does not exist: will be created"
    mkdir -p "${CACHEDIR}"
fi

#------------------------------------------------------------
# process each PDBID in turn
echo "processing PDBIDS: ${PDBIDS}"

# loop over all pdb ids
for pdbid in ${PDBIDS}; do
    # construct filename (lower case!)
    pdbfile=$(echo ${pdbid} | tr '[A-Z]' '[a-z]').pdb.gz
    url=${PDBURL}/${pdbfile}
    
    echo "[${pdbid}] PDB: ${pdbid} --> ${pdbfile}"

    # destination file
    pdbgz="${CACHEDIR}/${pdbfile}"
    pdb="${pdbgz%%.gz}"

    # check if the file already exists
    if [ -e "${pdbgz}" ] || [ -e "${pdb}" ]; then
	echo "[${pdbid}] ${pdbfile} already downloaded into ${CACHEDIR}"
	echo "------------------------------------------------------------"
	continue
    fi

    # we don't have the PDB yet so download it
    echo "[${pdbid}] PDB: url: ${url}"
    echo "[${pdbid}] downloading and uncompressing: ${pdbfile} --> ${pdb}" 
    # download and uncompress in one step
    curl ${url} | gunzip -c > ${pdb}
    errostatus=$?
    if [ ${errostatus} -eq 0 ]; then
	# successful
	echo "[${pdbid}] completed"
    else
        warn "Failed download of ${pdbid}"
	# clean up because gunzip leaves an empty file behind
	rm -f "${pdb}"
    fi
    echo "------------------------------------------------------------"
done

exit
