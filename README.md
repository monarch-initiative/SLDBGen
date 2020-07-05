# IDGKG-gen

This repository is generated for housing scripts, generated as part of the IDG Hackathon, for parsing, transforming and generating KGs from various sources.

The goal is to build scripts that is capable of parsing,
- Pharos/TCRD
- Monarch
- Datasets for Synthetic Lethality


The recommended interchange format: networkx JSON serialization, cytoscape JSON or TSV/CSV


## Setup 
The recommended setup is with a virtual environment

```bash
virtualenv py3
source py3/bin/activate
pip install -r requirements.txt 
```

and before each use of the scripts,
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