from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Patidar2020Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        """
        Persistent R-loops (RNA-DNA hybrids with a displaced single-stranded DNA) create DNA damage and lead to genomic
        instability. The 5'-3'-exoribonuclease 2 (XRN2) degrades RNA to resolve R-loops and promotes transcription
        termination.  XRN2 depletion compromised cell survival after additional knockdown of specific DNA repair
        proteins, including PARP1. These are shown in Figure 3B-D.
        . In the NHEJ group, 53BP1, DNA-PKcs, XRCC4 and Ligase4; in the HR group, Rad51, Brca1, Brca2, RPA and XRCC3,
         and in the BER group, Ape1, Fen1, Xrcc1, and Ligase3 were found to be synthetic lethal with XRN2 depletion
        """
        pmid = "32859985"
        super().__init__(fname=fname, pmid=pmid)

    def get_sli(self, geneA_sym, geneA_id, geneB_sym, geneB_id, sl=True):
        # Gene A should be eigther XRN2
        # These genes had activating mutations in the cell lines
        sli = SyntheticLethalInteraction(gene_A_symbol=geneA_sym,
                                         species_id="10090",
                                         gene_A_id=geneA_id,
                                         gene_B_symbol=geneB_sym,
                                         gene_B_id=geneB_id,
                                         gene_A_pert=SlConstants.SH_RNA,
                                         gene_B_pert=SlConstants.SI_RNA,
                                         effect_type=SlConstants.N_A,
                                         effect_size=SlConstants.N_A,
                                         cell_line=SlConstants.N_A,
                                         cellosaurus_id=SlConstants.N_A,
                                         cancer_type=SlConstants.N_A,
                                         ncit_id=SlConstants.N_A,
                                         assay=SlConstants.CELL_VIABILITY_ASSAY,
                                         pmid=self.pmid,
                                         SL=sl)
        return sli

    def parse(self):
        # SMARCA4-ARID2, SMARCA4-ACTB and SMARCC1-SMARCC2.
        xrn2 = "XRN2"
        xrn2id = self.get_ncbigene_curie(xrn2)
        # 53BP1=TP53BP1, DNA-PKcs=PRKDC, XRCC4 and Ligase4=LIG4
        nhej_sli = {'TP53BP1', 'PRKDC', 'XRCC4', 'LIG4'}
        # transcription termination factors (K-H, PSF, p54(nrb), and p15RS) did not reduce the survival of shXRN2 cells
        # (K-H=???, PSF=SFPQ, p54(nrb)=NONO, and p15RS=RPRD1A)
        non_sli = {'SFPQ', 'NONO', 'RPRD1A'}
        # HR group, Rad51, Brca1, Brca2, RPA and XRCC3
        hr_sli = {'RAD51', 'BRCA1', 'BRCA2', 'RPA1', 'XRCC3'}
        #  BER group, Ape1=APEX1, Fen1=Fen1, Xrcc1=XRCC1, and Ligase3=LIG3
        ber_sli = {'APEX1', 'FEN1', 'XRCC1', 'LIG3'}
        sli_genes = nhej_sli.union(hr_sli).union(ber_sli)
        sli_list = []
        for geneB in sli_genes:
            if geneB not in self.entrez_dict:
                raise ValueError("Could not find geneB:", geneB)
            geneBid = self.get_ncbigene_curie(geneB)
            sli = self.get_sli(geneA_sym=xrn2, geneA_id=xrn2id, geneB_sym=geneB, geneB_id=geneBid, sl=True)
            sli_list.append(sli)
        for geneB in non_sli:
            if geneB not in self.entrez_dict:
                raise ValueError("Could not find geneB:", geneB)
            geneBid = self.get_ncbigene_curie(geneB)
            sli = self.get_sli(geneA_sym=xrn2, geneA_id=xrn2id, geneB_sym=geneB, geneB_id=geneBid, sl=False)
            sli_list.append(sli)


        return sli_list
