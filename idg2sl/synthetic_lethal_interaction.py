from .parsers.sl_constants import SlConstants

class SyntheticLethalInteraction:
    """
    Instances of this class represent a single synthetic lethality (SL) interaction
    together with data about the experiment that was used to show the SL. A series
    of ingest scripts will transform the data of each of the ca. 30 published
    experiments to this common datastructure
    """

    def __init__(self,
                 gene_A_symbol=None,
                 gene_A_id=None,
                 gene_B_symbol=None,
                 gene_B_id=None,
                 gene_A_pert=None,
                 gene_B_pert=None,
                 effect_type=None,
                 effect_size=None,
                 species_id="9606",
                 cell_line="",
                 cellosaurus_id="",
                 cancer_type="",
                 ncit_id="",
                 assay=None,
                 pmid=None,
                 background_dependency_status=SlConstants.N_A,
                 background_dependency_gene_symbol=SlConstants.N_A,
                 background_dependency_gene_id=SlConstants.N_A,
                 SL=None):
        """

        """
        if gene_A_symbol is None:
            raise ValueError("Need to pass gene A")
        if gene_A_id is None:
            raise ValueError("Need to pass gene A id")
        if gene_B_symbol is None:
            raise ValueError("Need to pass gene B")
        if gene_B_id is None:
            raise ValueError("Need to pass gene B id")
        if gene_A_pert is None:
            raise ValueError("Need to pass gene A perturbation")
        if gene_B_pert is None:
            raise ValueError("Need to pass gene B perturbation")
        if assay is None:
            raise ValueError("Need to pass assay")
        if pmid is None:
            raise ValueError("Need to pass pmid")
        if SL is None:
            raise ValueError("Need to pass True or False for SL")
        self.gene_A_symbol = gene_A_symbol
        self.gene_A_id = gene_A_id
        self.gene_B_symbol = gene_B_symbol
        self.gene_B_id = gene_B_id
        self.gene_A_pert = gene_A_pert
        self.gene_B_pert = gene_B_pert
        self.assay = assay
        self.pmid = pmid
        # The cell line data is not obligatory. If it is not passed, set it to the empty string
        self.species_id = species_id
        self.cell_line = cell_line
        self.cellosaurus_id = cellosaurus_id
        self.cancer_type = cancer_type
        self.ncit_id = ncit_id
        if effect_type is None or effect_size is None:
            self.effect_type = ""
            self.effect_size = ""
        else:
            self.effect_type = effect_type
            self.effect_size = effect_size
        self._background_dependency_status = background_dependency_status
        self._background_dependency_gene_symbol = background_dependency_gene_symbol
        self._background_dependency_gene_id = background_dependency_gene_id
        self.SL = SL # True: synthetic lethal, False: negative control
        self.maximum_value = False

    def get_gene_A_symbol(self):
        return self.gene_A_symbol

    def get_gene_A_id(self):
        return self.gene_A_id

    def get_gene_B_symbol(self):
        return self.gene_B_symbol

    def get_gene_B_id(self):
        return self.gene_B_id

    def get_gene_A_pert(self):
        return self.gene_A_pert

    def get_gene_B_pert(self):
        return self.gene_B_pert

    def get_assay(self):
        return self.assay

    def get_pmid(self):
        return self.pmid

    def get_species_id(self):
        return self.species_id

    def get_cell_line(self):
        return self.cell_line

    def get_cellosaurus_id(self):
        return self.cellosaurus_id

    def get_cancer_type(self):
        return self.cancer_type

    def get_ncit_id(self):
        return self.ncit_id

    def get_effect_type(self):
        return self.effect_type

    def get_effect_size(self):
        return self.effect_size

    def set_maximum(self):
        self.maximum_value = True

    def is_maximum(self):
        return self.maximum_value

    def get_SL(self):
        return self.SL

    def get_tsv_line(self):
        if self.SL:
            sl = 'T'
        else:
            sl = 'F'
        lst = [self.gene_A_symbol,
               str(self.gene_A_id),
               self.gene_B_symbol,
               str(self.gene_B_id),
               self.gene_A_pert,
               self.gene_B_pert,
               self.effect_type,
               str(self.effect_size),
               self.species_id,
               self.assay,
               self.cell_line,
               self.cellosaurus_id,
               self.cancer_type,
               self.ncit_id,
               self.pmid,
               sl]
        return "\t".join(lst)

    @staticmethod
    def get_tsv_with_ensembl_header():
        lst = ['geneA',
               'geneA.ncbi-id',
               'geneA.ensembl-id',
               'geneB',
               'geneB.ncbi-id',
               'geneB.ensembl-id',
               'geneA.perturbation',
               'geneB.perturbation',
               'effect.type',
               'effect.size',
               'assay',
               'cell.line',
               'cellosaurus.id',
               'cancer',
               'ncit.id',
               'background.dependency.status',
               'background.dependency.genesymbol',
               'background.dependency.geneid',
               'pmid',
               'SL']
        return "\t".join(lst)

    def get_tsv_line_with_ensembl(self, symbol2ensembl):
        ensemblA = symbol2ensembl.get(self.gene_A_symbol, "n/a")
        ensemblB = symbol2ensembl.get(self.gene_B_symbol, "n/a")
        if self.gene_A_pert is None:
            print("NONE gene_A_pert")
            return ""
        if self.SL:
            sl = 'T'
        else:
            sl = 'F'
        lst = [self.gene_A_symbol,
               str(self.gene_A_id),
               ensemblA,
               self.gene_B_symbol,
               str(self.gene_B_id),
               ensemblB,
               self.gene_A_pert,
               self.gene_B_pert,
               self.effect_type,
               str(self.effect_size),
               self.assay,
               self.cell_line,
               self.cellosaurus_id,
               self.cancer_type,
               self.ncit_id,
               self._background_dependency_status,
               self._background_dependency_gene_symbol,
               self._background_dependency_gene_id,
               self.pmid,
               sl]
        return "\t".join(lst)


    def __repr__(self):
        """
        Returns:
            a string to be used by 'print'
        """
        return self.get_tsv_line()