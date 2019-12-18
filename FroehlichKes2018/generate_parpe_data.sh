#!/usr/bin/env bash
SCRIPT_PATH=$(dirname $BASH_SOURCE)
PARPE_DIR=/home/dweindl/src/parPE_2/

cd ${SCRIPT_PATH}

${PARPE_DIR}/misc/generateHDF5DataFileFromText.py \
    -o parpe_data.h5 \
    -s PEtab/CS_Signalling_ERBB_RAS_AKT_petab.xml \
    -d CS_Signalling_ERBB_RAS_AKT_petab \
    -m PEtab/measurements_petab.tsv \
    -c PEtab/conditions_petab.tsv \
    -n CS_Signalling_ERBB_RAS_AKT_petab \
    -p PEtab/parameters_petab.tsv
