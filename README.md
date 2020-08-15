# Synthetic Lethality Data Resource

This repository contains data derived from publications about
synthetic lethality in cancer together with Python code to
extract and transform the data into a common format and to
create a TSV output file.

To run the script, just enter
```buildoutcfg
$ python parse_human_SLI.py
```
The script will download the file ``protein-coding gene.txt``
from HGNC, which it uses to find NCBI Gene ids and Ensembl ids. 
We have extract relevant data from publications about synthetic
lethality (mainly from Supplemental Tables etc.). The script will 
create an output file called ``SL_data.tsv`` with positive and
negative (i.e., excluded) synthetic lethal interactions.


## Setup 
The package has a few requirements. The easiest way to set things
up is to use a virtual environment.

```bash
virtualenv py3
source py3/bin/activate
pip install -r requirements.txt 
```

and before each use of the script:
```bash
source py3/bin/activate
```

## Testing
Activate the virtual environment as above, and then install the nose package
```bash
source py3/bin/activate
pip install nose
nosetests
```