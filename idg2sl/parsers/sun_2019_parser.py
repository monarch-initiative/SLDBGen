from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Sun2019Parser(SL_DatasetParser):
    """
    A new sgRNA library was constructed for validation of the 926 candidate VHL SL genes by selecting ten new
    orthogonal sgRNAs per gene, based on an improved sgRNA picking algorithm (11,190 sgRNAs targeting the 926
    primary screening hits as well as 1,890 negative control sgRNAs)
    Together, these strategies allowed us to define a set of 350 VHL SL genes that were both expressed and
    functional in A-498 and 786-O cell lines (Fig. 2D and Supplemental Table S9).
    By applying fragments per kilobase of transcript per million mapped reads (FPKM) of less than 0.1 as a cutoff,
    36 (6%) and 14 (4%) genes for A-498 and 786-O cell lines, respectively, were eliminated (Supplemental Table S8).
    Together, these strategies allowed us to define a set of 350 VHL SL genes that were both expressed and functional
    in A-498 and 786-O cell lines
    For this script, I took the intersection of the genes listed in the three sheets of the original Supplemental Table S8
    and put them into a new file called 'data/Sun-TableS9.txt'
    """
    def __init__(self, fname='data/Sun-TableS9.txt'):
        pmid = '31436504'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        vhl = 'VHL'
        sli_list = []
        unclear_gene_symbols = {'QARS', 'SARS' }
        # I could figure out that the following mappings are correct and unique with the HGNC website
        mappings = {'ORAOV1': 'LTO1', 'VWA9': 'INTS14', 'NARFL':'CIAO3', 'WBSCR22': 'BUD23',
                    'UFD1L': 'UFD1', 'C7orf26': 'INTS15'}
        with open(self.fname) as f:
            for line in f:
                geneBsym = line.strip()
                if geneBsym in self.entrez_dict:
                    geneBid = self.get_ncbigene_curie(geneBsym)
                elif geneBsym == 'DARS' or geneBsym == 'NARS' or geneBsym == 'KARS' or geneBsym == 'YARS':
                    # A group of tRNA genes that need to have the '1' (I could map these uniquely with HGNC)
                    geneBsym = "%s1" % geneBsym
                    geneBid = self.get_ncbigene_curie(geneBsym)
                elif geneBsym in mappings:
                    geneBsym = mappings.get(geneBsym)
                    geneBid = self.get_ncbigene_curie(geneBsym)
                elif geneBsym in unclear_gene_symbols:
                    continue
                else:
                    raise ValueError("Could not find id for %s in Sun 2019" % geneBsym)
                sli = SyntheticLethalInteraction(gene_A_symbol=vhl,
                                                 gene_A_id=SlConstants.VHL_GENE_ID,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneBid,
                                                 gene_A_pert=SlConstants.LOF_MUTATION,
                                                 gene_B_pert=SlConstants.CRISPR_CAS9,
                                                 effect_type=SlConstants.N_A,
                                                 effect_size=SlConstants.N_A,
                                                 cell_line=SlConstants.A498_CELL,
                                                 cellosaurus_id=SlConstants.A498_CELLOSAURUS,
                                                 cancer_type=SlConstants.N_A,
                                                 ncit_id=SlConstants.N_A,
                                                 assay=SlConstants.CRISPR_CAS9_INTERFERENCE_ASSAY,
                                                 pmid=self.pmid,
                                                 SL=True)
                sli_list.append(sli)
        return sli_list
