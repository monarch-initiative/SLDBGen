from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Luo2009Parser(SL_DatasetParser):
    """
    Luo J, et al A genome-wide RNAi screen identifies multiple synthetic lethal
    interactions with the Ras oncogene. Cell. 2009 May 29;137(5):835-48.
    PubMed PMID: 19490893
    Supplemental Table S3 was copied from the original PDF file and stored to
    a file 'luo2009.tsv' in the data directory. The experiment used cell lines
    with activating (oncogenic) mutations in KRAS and used a first screen with
    relative abundance of shRNAs to identify 368 genes using stringent criteria.
    They tested 320 candidates from the first screen and found 83 shRNAs targeting
    77 genes to preferentially decreased the viability of the KRAS mutant cells
    compared to WT cells, thus indicating SL. They screened 68 of these candidates
    in a second line and could confirm 50 of them (73.5%), indicating that the
    majority of the candidates were likely to be true positives.
    Note -- In our data file, I changed QARS to QARS1, which is the current symbol of NM_005051
    I also changed UCRC to UQCR10, which is the current symbol of NM_001003684
    """

    def __init__(self, fname='data/luo2009.tsv'):
        """
        data/luo2009.tsv is our copy of Supplemental Table S3
        """
        pmid = '19490893'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        kras_symbol = 'KRAS'
        kras_id = 'NCBIGene:3845'
        kras_perturbation = SlConstants.ACTIVATING_MUTATION
        gene2_perturbation = SlConstants.SH_RNA
        assays = ['competitive hybridization', 'multicolor competition assay']
        assay_string = ";".join(assays)
        effect_type = 'stddev'
        cell_line = SlConstants.DLD1_CELL
        cellosaurus = SlConstants.DLD1_CELLOSAURUS
        cancer = SlConstants.COLORECTAL_CARCINOMA
        ncit = SlConstants.COLORECTAL_CARCINOMA_NCIT

        # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
        # Symbol	Accession	v2SH	Sequence	Mean.DLD1	SD.DLD1	Mean.HCT116	SD.HCT116	
        sli_dict = defaultdict(list)
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) < 8:
                    raise ValueError("Bad line in Luo2009 with less than 8 fields")
                geneB_sym = self.get_current_symbol(row['Symbol'])
                if geneB_sym == 'CXORF40A':
                    geneB_sym = 'EOLA1'
                if geneB_sym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneB_sym))
                elif geneB_sym == 'FLJ34747':
                    # This is a LINC, plus the symbol is old
                    geneB_sym = 'LINC00547'
                    geneB_id = 'NCBIGene:400121'
                elif geneB_sym == 'LOC283194':
                    geneB_id = 'NCBIGene:283194' # an ncRNA
                elif geneB_sym == 'LOC285556':
                    geneB_id = 'NCBIGene:285556'
                elif geneB_sym == 'LOC149654' or geneB_sym == 'LOC730000':
                    continue  # Could not find these in NCBI Gene or HCNG
                else:
                    raise ValueError("Could not get NCBI id for gene %s in Luo2009" % geneB_sym)
                stddev = float(row['SD.DLD1'])  # float(fields[5])
                SL = True  # All data in this set is True
                sli = SyntheticLethalInteraction(gene_A_symbol=kras_symbol,
                                                 gene_A_id=kras_id,
                                                 gene_B_symbol=geneB_sym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=kras_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type=effect_type,
                                                 effect_size=stddev,
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay_string,
                                                 pmid=self.pmid,
                                                 SL=SL)
                gene_pair = GenePair(kras_symbol, geneB_sym)
                sli_dict[gene_pair].append(sli)
        sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
