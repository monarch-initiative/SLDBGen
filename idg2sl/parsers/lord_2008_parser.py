from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Lord2008Parser(SL_DatasetParser):
    def __init__(self, fname='data/lord-PARP1-2008.tsv'):
        pmid = 'PMID:18832051'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        """
        Parsing data from
        Lord CJ,et al A high-throughput RNA interference screen for DNA repair determinants of PARP inhibitor sensitivity.
        DNA Repair (Amst). 2008 Dec 1;7(12):2010-9. PMID: 18832051.
        This siRNA library targeted over 98% of all of the known DNA repair proteins, as defined by Wood et al. [16] and
        in particular encompassed siRNA targeting all the major components of BER (21 genes), MMR (11 genes),
        NER (28 genes), HR(19 genes) and NHE
        The HTS was performed using the diploid CAL51 breast cancer cell line.
        We used a log 2 surviving fraction (log 2 SF) threshold of −0.1 or less to identify siRNAs that significantly
        sensitized to KU0058948 in the HTS (Supplementary Tables 1–3). This threshold represented the limit of three
        standard deviations from the median of effects in siCON transfected cells. 67 of the 460 experimental siRNA
        satisfied this hit criteria (Fig. 3A).
        Using this distinction, eight genes, plus the control (BRCA1) were identified where both siRNA in the library.
        significantly sensitized to KU0058948 (Fig. 3C); ATR, BRCA2, DDB1, LIG1, PCNA, RAD51, XAB2 and XRCC1.
        We have previously reported, using a different assay system, that silencing of ATR, BRCA2 or RAD51 expression
        significantly sensitizes cells to KU0058948, presumably by causing defective HR [3,6].
        The following parse takes these 8 interactions to be true SL. We additionally generate a list of negative SLs
        using the criteria that neither of the cells attained an effect that was even half as strong
        # The data come as pairs of lines. We will first store the pairs and then only take those that satisfy the
        criteria for positive or negative (PARP sensitization). There are actually 9 SLA pairs including RPA3, which
        is shown in Figure 4C but not mentioned in the text. That is what we get with this parse!
        """
        parp1_symbol = 'PARP1'
        parp1_id = 'NCBIGene:142'
        parp1_perturbation = SlConstants.PHARMACEUTICAL.to_string()
        gene2_perturbation = SlConstants.SI_RNA.to_string()

        assay_string = SlConstants.RNA_INTERFERENCE_ASSAY.to_string()
        effect_type = 'stddev'
        cell_line = 'CAL-51'
        cellosaurus = 'CVCL_1110'
        cancer = "Breast Carcinoma"
        ncit = "NCIT:C4872"
        parpdict = defaultdict(list)
        # The following list includes symbols that are not current but either could
        # not be matched or match to multiple possible candidates (PMS2L4 is a pseudogene)
        unclear_gene_symbols = {'CDC2', 'NBS1', 'TGIF', 'PMS2L4'}
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) != 3:
                    raise ValueError("Malformed line with length %d instead of 3: %s" % (len(row), row))
                geneBsym = row['gene']
                if geneBsym in unclear_gene_symbols:
                    continue
                geneBsym = self.get_current_symbol(geneBsym)
                parp_sens = float(row['parp_sens'])
                if geneBsym == 'BRCA1':
                    continue
                elif geneBsym == 'GFP-22' or geneBsym == 'SCRAM':
                    continue  # a control siRNA
                # ignore the third field
                parpdict[geneBsym].append(parp_sens)
        sli_list = []
        for geneBsym, parp_sens_list in parpdict.items():
            if geneBsym in self.entrez_dict:
                geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneBsym))
            else:
                raise ValueError("Could not find iid for %s in Lord 2008 " % geneBsym)
            if len(parp_sens_list) != 2:
                raise ValueError("Length of list not equal to 2 for %s (len was %d)" % (
                    geneBsym, len(parp_sens_list)))  # should never happen
            if parp_sens_list[0] <= -0.1 and parp_sens_list[1] <= -0.1:
                sli = SyntheticLethalInteraction(gene_A_symbol=parp1_symbol,
                                                 gene_A_id=parp1_id,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=parp1_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type=effect_type,
                                                 effect_size=min(parp_sens_list[0], parp_sens_list[1]),
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay_string,
                                                 pmid=self.pmid,
                                                 SL=True)
                sli_list.append(sli)
            else:
                effectsize = min(parp_sens_list[0], parp_sens_list[1])
                if effectsize > -0.05:
                    sli = SyntheticLethalInteraction(gene_A_symbol=parp1_symbol,
                                                     gene_A_id=parp1_id,
                                                     gene_B_symbol=geneBsym,
                                                     gene_B_id=geneB_id,
                                                     gene_A_pert=parp1_perturbation,
                                                     gene_B_pert=gene2_perturbation,
                                                     effect_type=effect_type,
                                                     effect_size=min(parp_sens_list[0], parp_sens_list[1]),
                                                     cell_line=cell_line,
                                                     cellosaurus_id=cellosaurus,
                                                     cancer_type=cancer,
                                                     ncit_id=ncit,
                                                     assay=assay_string,
                                                     pmid=self.pmid,
                                                     SL=False)
                    sli_list.append(sli)
        return sli_list