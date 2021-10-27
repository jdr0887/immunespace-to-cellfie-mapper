# immunespace-to-cellfie-mapper


To generate the genename-to-hgnc-id.csv file:
```shell
curl http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt | \
  awk '{print $2 "," $1}' | sed 's;HGNC:;;g' > genename-to-hgnc-id.csv && sed -i 1d genename-to-hgnc-id.csv
```
