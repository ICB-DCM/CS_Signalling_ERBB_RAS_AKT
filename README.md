# CS_Signalling_ERBB_RAS_AKT

Model and data from FroehlichKes2018
https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30438-1 

and 

model and data from SchmiesterSch2019 
https://www.biorxiv.org/content/10.1101/579045v1


## Contents

- `generate_amici_model.sh`

  Generate AMICI model CS_Signalling_ERBB_RAS_AKT from PEtab files 


- `generate_parpe_data.sh`

  Generate parPE training data HDF5 file for CS_Signalling_ERBB_RAS_AKT from 
  PEtab files. Requires output of `generate_amici_model.sh`.

- `hierarchical-optimization-manuscript/`

  Model and data from SchmiesterSch2019 
  https://www.biorxiv.org/content/10.1101/579045v1

- `legacy-format`

   CS_Signalling_ERBB_RAS_AKT in legacy format with import scripts to generate
   PEtab-format model and data as in `PEtab/`
