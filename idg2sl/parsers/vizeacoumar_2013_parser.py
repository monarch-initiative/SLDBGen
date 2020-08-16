from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Vizeacoumar2013Parser(SL_DatasetParser):
    """
    six isogenic cell lines were screened  in parallel using a standardized genome-scale pooled shRNA screening pipeline.
    The HCT116 genetic background was chosen because it is near diploid with intact DNA damage and spindle checkpoints.
    The ‘query' genotypes chosen were PTTG1−/−, BLM−/−, MUS81−/−, PTEN−/− and KRAS+/−
    """

    def __init__(self, fname=None):
        pmid = 'PMID:24104479'
        super().__init__(fname=fname, pmid=pmid)
        self.sli_list = []
        # old symbols that are assigned to multiple genes
        # Pseudogenes: PMS2L5
        self.unclear_gene_symbols = {'51639', 'PMS2L5'}

    def parseLoF(self, geneA, fname):
        """
        BLM, MUS81, PTEN, PTTG1
        """
        geneAid = self.entrez_dict.get(geneA)
        geneA_perturbation = SlConstants.LOF_MUTATION
        gene2_perturbation = SlConstants.SI_RNA
        assay_string = SlConstants.MULTICOLOR_COMPETITION_ASSAY
        cell_line = 'HCT 116'
        cellosaurus = 'CVCL_0291'
        cancer = SlConstants.N_A
        ncit = SlConstants.N_A
        c = 0
        with open(fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if not row['Expression'] == 'Expressed':
                    continue
                geneBsym = self.get_current_symbol(row['human gene'])
                if geneBsym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneBsym))
                elif geneBsym in self.unclear_gene_symbols:
                    continue
                else:
                    raise ValueError("Could not find iid for %s in Brough 2018 2008 " % geneBsym)
                conf80 = int(row['80% Confidence Interval (P<0.2)'])
                if conf80 == 1:
                    c += 1
                    sli = SyntheticLethalInteraction(gene_A_symbol=geneA,
                                                 gene_A_id=geneAid,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=geneA_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type='confidence.80%',
                                                 effect_size='true',
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay_string,
                                                 pmid=self.pmid,
                                                 SL=True)
                    self.sli_list.append(sli)
        print("%s got %d genes" % (geneA, c))

    def parseKRAS(self):
        """
        BLM is BLM RecQ like helicase
        """
        geneA = 'KRAS'
        geneAid = self.get_current_symbol(geneA)
        fname = 'data/vizeacoumarSuppl4-PTEN.tsv'
        geneA_perturbation = SlConstants.ACTIVATING_MUTATION
        gene2_perturbation = SlConstants.SI_RNA
        assay_string = SlConstants.MULTICOLOR_COMPETITION_ASSAY
        cell_line = 'HCT 116'
        cellosaurus = 'CVCL_0291'
        cancer = SlConstants.N_A
        ncit = SlConstants.N_A
        c = 0
        with open(fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if not row['Expression'] == 'Expressed':
                    continue
                geneBsym = self.get_current_symbol(row['human gene'])
                if geneBsym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneBsym))
                elif geneBsym in self.unclear_gene_symbols:
                    continue
                else:
                    raise ValueError("Could not find iid for %s in Brough 2018 2008 " % geneBsym)
                conf80 = int(row['80% Confidence Interval (P<0.2)'])
                if conf80 == 1:
                    c += 1
                    sli = SyntheticLethalInteraction(gene_A_symbol=geneA,
                                                 gene_A_id=geneAid,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=geneA_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type='confidence.80%',
                                                 effect_size='true',
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay_string,
                                                 pmid=self.pmid,
                                                 SL=True)
                    self.sli_list.append(sli)
        print("%s got %d genes" % (geneA, c))

    def parse(self):
        blm = 'BLM'
        blm_fname = 'data/vizeacoumarSuppl4-BLM.tsv'
        self.parseLoF(geneA=blm, fname=blm_fname)
        mus81 = 'MUS81'
        mus81_fname = 'data/vizeacoumarSuppl4-MUS81.tsv'
        self.parseLoF(geneA=mus81, fname=mus81_fname)
        pttg1 = 'PTTG1'
        pttg1_fname = 'data/vizeacoumarSuppl4-PTTG1.tsv'
        self.parseLoF(geneA=pttg1,  fname=pttg1_fname)
        pten = 'PTEN'
        pten_fname = 'data/vizeacoumarSuppl4-PTEN.tsv'
        self.parseLoF(geneA=pten,  fname=pten_fname)
        self.parseKRAS()
        return self.sli_list

