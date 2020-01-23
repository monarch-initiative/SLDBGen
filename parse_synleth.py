import os, re, sys
import pandas as pd
import networkx as nx
import json

human_gene_lookup = None
def map_symbol_to_gene(symbol):
	global human_gene_lookup
	if human_gene_lookup is None:
		human_gene_lookup = {}
		data = json.load(open('lookup/monarch_human_gene_lookup.json'))
		docs = data['response']['docs']
		for d in docs:
			human_gene_lookup[d['subject_label']] = d['subject']

	retval = symbol
	if symbol in human_gene_lookup:
		retval = human_gene_lookup[symbol]
	return retval


table = pd.read_csv(sys.argv[1], sep='\t')

graph = nx.MultiDiGraph()

for index, row in table.iterrows():
	n1 = row['gene1']
	n1_curie = map_symbol_to_gene(n1)
	n2 = row['gene2']
	n2_curie = map_symbol_to_gene(n2)
	n1_props = {'gene1 perturbation': row['gene1 perturbation']}
	n2_props = {'gene2 perturbation': row['gene2 perturbation']}
	publication = f'PMID:{row["PMID"]}'
	edge_props = {'edge_label': 'regulates', 'relation': 'RO:0002211', 'publications': [publication], 'cancer_type_tested': row['cancer type tested'], 'synthetic_lethality': bool(row['SL'])}
	graph.add_node(n1_curie, **n1_props)
	graph.add_node(n2_curie, **n2_props)
	graph.add_edge(n1_curie, n2_curie, **edge_props)

OUTSTREAM = open('output.cytoscape.json', 'w')
nx_data = nx.cytoscape_data(graph)
nx_data_json = json.dumps(nx_data)
OUTSTREAM.write(nx_data_json)
