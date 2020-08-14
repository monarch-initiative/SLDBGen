from collections import defaultdict
from idg2sl.gene_pair import GenePair
from idg2sl.synthetic_lethal_interaction import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
import csv


class Turner2008Parser(SL_DatasetParser):
    """
    Parse data from  Turner NC, et al., A synthetic lethal siRNA screen identifying genes mediating
    sensitivity to a PARP inhibitor. EMBO J. 2008 May 7;27(9):1368-77. PubMed PMID: 18388863
    PARP1 was inhibited by the PARP inhibitor KU0058948, and a short interfering RNA library targeting
    779 human protein kinase and kinase assocaited genes was applied.
    :param path:
    I replaced MGC5601 by ACAD10 (Gene:80724)
    HSMDPKIN by CDC42BPG Gene ID: 55561,
    FLJ35107 by TPRXL Gene ID: 348825,
    """

    def __init__(self, fname='data/turner-PARP1-2008.tsv'):
        pmid = 'PMID:18388863'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        parp1_symbol = 'PARP1'
        parp1_id = 'NCBIGene:142'
        parp1_perturbation = SlConstants.DRUG
        gene2_perturbation = SlConstants.SI_RNA
        assays = ['competitive hybridization', 'multicolor competition assay']
        assay_string = ";".join(assays)
        effect_type = 'stddev'
        cell_line = 'CAL-51'
        cellosaurus = 'CVCL_1110'
        cancer = "Breast Carcinoma"
        ncit = "NCIT:C4872"
        sli_dict = defaultdict(list)
        with open(self.fname) as csvfile:
            # SMARTpool	Z score	percent-siCONTROL
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) < 3:
                    raise ValueError("Bad row for Turner et al, only %d fields found (%s)" % (len(row), row))
                geneB_sym = self.get_current_symbol(row['SMARTpool'])
                zscore = float(row['Z score'])
                if geneB_sym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneB_sym))
                elif geneB_sym == 'IMPK':
                    continue  # could not be found in HGNC or NCBI Gene
                elif geneB_sym == 'FLJ34389':
                    geneB_sym = 'MLKL'
                    geneB_id = "NCBIGene:197259"
                else:
                    if zscore > 3.0:
                        raise ValueError("Could not get NCBI id for gene %s in Turner 2008" % geneB_sym)
                    else:
                        continue # These are negative examples, we will just skiip
                if zscore <= -3.0:
                    SL = True
                else:
                    SL = False
                sli = SyntheticLethalInteraction(gene_A_symbol=parp1_symbol,
                                                 gene_A_id=parp1_id,
                                                 gene_B_symbol=geneB_sym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=parp1_perturbation.to_string(),
                                                 gene_B_pert=gene2_perturbation.to_string(),
                                                 effect_type=effect_type,
                                                 effect_size=zscore,
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay_string,
                                                 pmid=self.pmid,
                                                 SL=SL)
                gene_pair = GenePair(parp1_symbol, geneB_sym)
                sli_dict[gene_pair].append(sli)
                sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
