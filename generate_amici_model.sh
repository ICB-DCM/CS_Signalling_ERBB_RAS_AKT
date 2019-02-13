#!/usr/bin/env bash
# Generate AMICI model
amici_import_petab.py -v -s 'PEtab/CS_Signalling_ERBB_RAS_AKT_petab.xml' \
    -c 'PEtab/conditions_petab.tsv'
    -m 'PEtab/measurements_petab.tsv'
    -p 'PEtab/parameters_petab.tsv'

