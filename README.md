# *C*ontextual *H*ub *A*nalysis *T*ool
CHAT identifies nodes in a network that are significantly more connected to contextually relevant nodes (e.g. differentially expressed genes) than is expected by chance.

###To run CHAT:

Provide an XGMML formatted network file (see [chat.primesdb.eu](http://chat.primesdb.eu/cli.html) for format requirements) and specify the attribute containing context:
```
chat.py -x filename -a attributeName
```
OR

A tab delimited text file of human/mouse Ensembl Gene IDs (first column, labeled with "ens") plus contextual information (second column, labeled "pVal"). CHAT will build a network and then run the analysis:

```
chat.py -n filename
```

###Other arguments:
```
  -x xgmml              Input xgmml file
  -n genes              List of genes and context data for network
  -h                    Show this help message and exit
  -a attName            Attribute that contains context (case sensitive)
  -e ensAtt             Attribute that contains ENS ID (case sensitive)
  -pV pThresh           P-value threshold, default 0.1
  -fc foldChangeBounds  Fold-change boundaries, e.g. 2
  -g hubDeGree          Number of connections to classify a hub (default n=5)
  -m mitab              MITAB file of interactions
  -p prefix             Specify a prefix to be appended to output files.
  -v                    Print all hub parameters
  ```
