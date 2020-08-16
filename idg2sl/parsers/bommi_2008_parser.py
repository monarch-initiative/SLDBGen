from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Bommi2008Parser(SL_DatasetParser):
    """
    Bommi-Reddy A, et al  Kinase requirements in human cells: III. Altered kinase
    requirements in VHL-/- cancer cells detected in a pilot synthetic lethal screen.
    Proc Natl Acad Sci U S A. 2008 Oct 28;105(43):16484-9. PubMed PMID: 18948595.
    The paper describes the use of two cell types: 786-O (Table S4) and RCC4 (Table S5).
    Each of these cells has a loss of function mutation in the tumor suppressor gene VHL.
    The defects in the other gene were induced with a lentiviral vector that produced a shRNA. The
    final results are in Tables S4 and S5.
    The authors require a differential viability above 35% for the 786-O cells or a differential
    viability above 20% for the RCC4 cells. We have copied that values from the supplemental tables
    (which were PDF files) into a TSV file called data/bommi-reddy-2008.tsv. For the output file,
    we demand that at least one shRNA was above threshold, and just report the best outcome per gene.
     The authors showed experimentally that two genes, IRR and HER4, actually are not SL, so we remove these
    by hand
    I replaced JNK3 by MAPK10
    I replaced JNK2 by MAPK9
    I replaced AMPKa1 by PRKAA1
    I prepleaced KHS1 by MAP4K5
    I replaced SgK495 by STK40
    I replaced SURTK106 by STYK1
    """

    def __init__(self, fname='data/bommi-reddy-2008.tsv'):
        pmid = 'PMID:18948595'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        # because of the experiment, geneA is always VHL.
        vhl_symbol = 'VHL'
        vhl_id = SlConstants.VHL_GENE_ID
        vhl_perturbation = SlConstants.LOF_MUTATION
        gene2_perturbation = SlConstants.SH_RNA
        assays = [SlConstants.COMPETITIVE_HYBRIDIZATION, SlConstants.MULTICOLOR_COMPETITION_ASSAY]
        assay_string = ";".join(assays)
        effect_type = 'differential_viability'
        cell_786O = "786-0"
        cellosaurus_786O = "CVCL_1051"
        cell_RCC4 = "RCC4"
        cellosaurus_RCC4 = "CVCL_0498"
        # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
        sli_dict = defaultdict(list)
        # The following list includes symbols that are not current but either could
        # not be matched or match to multiple possible candidates
        unclear_gene_symbols = {'PITSLRE', 'TAK1', 'PKD3', 'CAMLCK', 'MAPAPK3', 'CK1E', 'CK2A2', 'PDGRFB', 'ZC1/HGK'}
        # gene	differential	cell	table
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) < 4:
                    raise ValueError("Only got %d fields but was expecting 4" % len(row))
                genesy = row['gene'].upper()
                geneB_sym = self.get_current_symbol(genesy)
                if geneB_sym == "IRR" or geneB_sym == "HER4":
                    continue
                if geneB_sym in unclear_gene_symbols:
                    continue  # Symbol could be either CDK11A or CDK11B
                if geneB_sym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneB_sym))
                else:
                    raise ValueError("Could not find id for %s in Bommi 2008" % geneB_sym)
                effect = float(row['differential'])
                cell = row['cell']
                if cell == 'RCC4':
                    cell_line = cell_RCC4
                    cellosaurus = cellosaurus_RCC4
                elif cell == '786-0':
                    cell_line = cell_786O
                    cellosaurus = cellosaurus_786O
                else:
                    raise ValueError("Did not recognize cell type '%s'" % cell)
                table = row['table']
                assay_string = "differential viability assay {}({})".format(cell, table)
                SL = True  # All data in this set is True # TODO CHECK
                sli = SyntheticLethalInteraction(gene_A_symbol=vhl_symbol,
                                                 gene_A_id=vhl_id,
                                                 gene_B_symbol=geneB_sym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=vhl_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type=effect_type,
                                                 effect_size=effect,
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=SlConstants.CLEAR_CELL_RENAL_CELL_CARCINOMA,
                                                 ncit_id=SlConstants.CLEAR_CELL_RENAL_CELL_CARCINOMA_NCIT,
                                                 assay=assay_string,
                                                 pmid=self.pmid,
                                                 SL=SL)
                gene_pair = GenePair(vhl_symbol, geneB_sym)
                sli_dict[gene_pair].append(sli)
        sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
