import uuid

from biolink_model_pydantic.model import ( #type: ignore
    Gene,
    GeneToGeneAssociation,
    Predicate
)

from koza.cli_runner import koza_app #type: ignore

source_name="sl_data_config"

row = koza_app.get_row(source_name)

# Entities

g1 = Gene(id="ensembl:"+str(row['geneA.ensembl-id']),
            name=row["geneA"])
g2 = Gene(id="ensembl:"+str(row['geneB.ensembl-id']),
            name=row["geneB"])

# Association

association = GeneToGeneAssociation(
        id="uuid:" + str(uuid.uuid1()),
        subject=g1.id,
        predicate="biolink:genetically_interacts_with",
        object=g2.id,
    )

koza_app.write(g1, association, g2)