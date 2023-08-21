### AGORA2.01 - Download

Downloading all AGORA2 models in SBML format:

```sh
while read p; do
  echo "$p"
  wget https://www.vmh.life/files/reconstructions/AGORA2/version2.01/sbml_files/individual_reconstructions/$p -O sbml/$p
  gzip sbml/$p
done <sbml_files.txt
```

### SBML to 'modelorg' object in R

R-objects are stored in RDS files.

see script `xml2rds.R`

> A number of models from AGORA2 (v2.01) seem to be in an invalid SBML file format. In those cases, R-objects could not be created. The exact error messages from libSBML when trying to read the respective files is stored in the log file `errlog.txt`.

### Predict auxotrophies

Please see script `predict_auxotrophies_AGORA2.R`.

Auxotrophy predictions for each model are written in file `AGORA2_predicted_auxotrophies.csv`.
