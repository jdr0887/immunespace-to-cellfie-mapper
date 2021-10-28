# immunespace-to-cellfie-mapper


To generate the genename-to-hgnc-id.csv file:
```shell
curl http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt | \
  awk -F'\t' '{print $2 "," $1 "," $19}' | sed 's;HGNC:;;g' > genename-data.csv && sed -i 1d genename-data.csv
```
