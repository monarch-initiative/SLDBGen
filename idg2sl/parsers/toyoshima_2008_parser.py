from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Toyoshima2008Parser(SL_DatasetParser):
    def __init__(self, fname='data/toyoshima-MYC-2008.tsv'):
        """
        Parsing data from
        Toyoshima M, et al. Functional genomics identifies therapeutic targets for MYC-driven cancer.
        Proc Natl Acad Sci U S A. 2012 Jun 12;109(24):9545-50. PMID: 22623531
        The results of the screen revealed 148 hits, defined according to a Z score of ≥2 (23), including 140 genes and
        eight microRNAs (Fig. 1B). Here, we focus on the 140 gene hits, which we designate MYC-synthetic lethal (MYC-SL)
        genes. To eliminate siRNAs that exhibited substantial growth inhibition properties in normal cells, siRNAs
        with >50% reduced viability in HFF-pBabe were eliminated from further consideration regardless of differential
        toxicity. This process left 102 MYC-SL gene hits for follow-up
        NOTE: I see only 101 genes in the Supplemental table.
        """
        pmid = '22623531'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        mycsymbol = 'MYC'
        effect_type = 'stddev'
        cell_line = 'HFF-Myc'
        cellosaurus = 'CVCL_Y511'
        sl_list = []
        # The following list includes symbols that are not current but either could
        # not be matched or match to multiple possible candidates
        unclear_gene_symbols = {'MLCK'}
        # Gene.Symbol	Accession.number	Z.score.greaterthan	%Viability.HFF-pB	%Viability.HFF-MYC	Ratio pBabe/Myc
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) != 6:
                    raise ValueError("Bad line with %d fields: %s" % (len(row), row))
                # Gene Symbol	Accession number	Z score (>than)	%Viability HFF-pB	%Viability HFF-MYC	Ratio pBabe/Myc
                geneBsym = row['Gene.Symbol']
                if geneBsym in unclear_gene_symbols:
                    continue
                geneBsym = self.get_current_symbol(geneBsym)
                if geneBsym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneBsym))
                else:
                    raise ValueError("Could not find id for symbol %s in Toyoshima 2008" % geneBsym)
                zscore = float(row['Z.score.greaterthan'])
                sli = SyntheticLethalInteraction(gene_A_symbol=mycsymbol,
                                                 gene_A_id=SlConstants.MYC_GENE_ID,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=SlConstants.OVEREXPRESSION,
                                                 gene_B_pert=SlConstants.SI_RNA,
                                                 effect_type=effect_type,
                                                 effect_size=zscore,
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=SlConstants.N_A,
                                                 ncit_id=SlConstants.N_A,
                                                 assay=SlConstants.RNA_INTERFERENCE_ASSAY,
                                                 pmid=self.pmid,
                                                 SL=True)
                sl_list.append(sli)
        return sl_list
