import uuid

from biolink_model_pydantic.model import ( #type: ignore
    Gene,
    GeneToGeneAssociation
)

from koza.cli_runner import koza_app #type: ignore

source_name="sl_data_config"

# CURIEs for perturbation types
perturb_lookup = {"CRISPR CAS9":"BAO:0010249",
                    "RNA-interference assay":"BAO:0002434",
                    "activating mutation":"BAO:0000087",
                    "agonist":"BAO:0002902",
                    "antisense oligonucleotide":"BAO:0002439",
                    "cohort study":"SIO:001067",
                    "degradation":"INO:0000076",
                    "inhibitory antibody":"BAO:0000091",
                    "knockout":"BAO:0002441",
                    "lof_mutation":"NCIT:C178119",
                    "natural (is a TSG)":"NCIT:C36334",
                    "overexpression":"BAO:0002442",
                    "pharmaceutical":"CHEBI:52217",
                    "promoter hypermethylation":"NCIT:C20102",
                    "sgRNA":"SO:0001998",
                    "shRNA":"BAO:0000323",
                    "siRNA":"BAO:0000324",}

full_assay_lookup = {"Big Papi":"",
                    "CRISPR-Cas9 Interference assay":"BAO:0010249",
                    "PDX":"",
                    "RNA-interference assay":"BAO:0002434",
                    "apoptosis assay":"BAO:0002043",
                    "cell viability assay":"BAO:0003009",
                    "cisplatin toxicity assay":"BAO:0000623",
                    "competitive hybridization;multicolor competition assay":"",
                    "cytotoxicity assay":"",
                    "differential viability assay":"",
                    "growth inhibition assay":"BAO:0000094",
                    "multicolor competition assay":"",
                    "pharmaceutical + siRNA":"",
                    "pharmaceutical inhibition assay":"",
                    "proportions.of.sense.and.antisense.insertions":"",
                    "sgRNA depletion assay":"",
                    "shRNA":"BAO:0000323",
                    "shRNA depletion assay":"BAO:0000323",
                    "siMEM+penetrance":"",
                    "transgenic mouse model":"NCIT:C14240"}

row = koza_app.get_row(source_name)

# Entities

g1 = Gene(id="ENSEMBL:"+str(row['geneA.ensembl-id']),
            name=row["geneA"],
            xref=[row['geneA.ncbi-id']])
g2 = Gene(id="ENSEMBL:"+str(row['geneB.ensembl-id']),
            name=row["geneB"],
            xref=[row['geneB.ncbi-id']])

# Association

association = GeneToGeneAssociation(
        id="uuid:" + str(uuid.uuid1()),
        subject=g1.id,
        predicate="biolink:genetically_interacts_with",
        object=g2.id,
        publications=["PMID:"+str(row['pmid'])],
        qualifiers=[]
    )
if perturb_lookup[str(row["geneA.perturbation"])] != "":
    association.qualifiers.append(perturb_lookup[str(row["geneA.perturbation"])])
else:
    association.qualifiers.append(str(row["geneA.perturbation"]))
if perturb_lookup[str(row["geneB.perturbation"])] != "":
    association.qualifiers.append(perturb_lookup[str(row["geneB.perturbation"])])
else:
    association.qualifiers.append(str(row["geneB.perturbation"]))

if full_assay_lookup[str(row["assay"])] != "":
    association.qualifiers.append(full_assay_lookup[str(row["assay"])])
else:
    association.qualifiers.append(str(row["assay"]))

if str(row["cell.line"]) != "n/a":
    if str(row["cellosaurus.id"]) != "n/a":
        association.qualifiers.append("cellosaurus:"+str(row["cellosaurus.id"]))
    else:
        association.qualifiers.append(str(row["cell.line"]))

koza_app.write(g1, association, g2)