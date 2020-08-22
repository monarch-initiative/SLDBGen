from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Martin2010and2011Parser(SL_DatasetParser):
    """
    This class ingests data from PMID: 20227038 and PMID: 21242281, both by the same group
    """
    def __init__(self, fname=None):
        pmid = "20227038"
        super().__init__(fname=fname, pmid=pmid)

    def create_sli(self, geneA, geneB, cell, cellosaurus, pmid, SL=True):
        """
        MLH1/POLG and MSH2/POLB SSLs Are Rescued by MUTYH Silencing
        """
        geneAid = self.get_ncbigene_curie(geneA)
        geneBid = self.get_ncbigene_curie(geneB)
        mutyh = 'MUTYH'
        mutyh_id = self.get_ncbigene_curie(mutyh)
        sli = SyntheticLethalInteraction(gene_A_symbol=geneA, gene_A_id=geneAid,
                                         gene_B_symbol=geneB, gene_B_id=geneBid,
                                         gene_A_pert=SlConstants.LOF_MUTATION,
                                         gene_B_pert=SlConstants.SI_RNA,
                                         cell_line=cell,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=SlConstants.N_A,
                                         ncit_id=SlConstants.N_A,
                                         effect_size=SlConstants.N_A, effect_type=SlConstants.N_A,
                                         assay=SlConstants.CELL_VIABILITY_ASSAY,
                                         background_dependency_status=SlConstants.WILDTYPE,
                                         background_dependency_gene_symbol=mutyh,
                                         background_dependency_gene_id=mutyh_id,
                                         SL=SL,
                                         pmid=pmid)
        return sli




    def martin_2011(self):
        """
        C9ORF96 is now STKLD1 (serine/threonine kinase like domain containing 1)
        DKFZP434C131 is now ULK3 (unc-51 like kinase 3)
        ICK is now CILK1 (ciliogenesis associated kinase 1)
        TRAP1 is TNF receptor-associated protein 1  (see PMID: 23525905)
        Parkin is PRKN
        FLJ23356 is now POMK (protein O-mannose kinase)
        """
        sli_list = []
        pmid = '21242281'
        mlh1 = 'MLH1'
        mlh1_sl_parters = {'STKLD1', 'ULK3', 'GUCY2D', 'CILK1', 'MAP2K6', 'ROR2', 'PINK1'}
        for geneB in mlh1_sl_parters:
            sli = self.create_sli(geneA=mlh1, geneB=geneB,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS, pmid=pmid)
            sli_list.append(sli)
        msh2_sl_partners = {'CKMT2', 'POMK', 'NEK11', 'PCK2', 'RPS6KB1', 'LTK', 'PINK1'}
        msh2 = 'MSH2'
        for geneB in msh2_sl_partners:
            sli = self.create_sli(geneA=msh2, geneB=geneB,
                              cell=SlConstants.HEC59_CELL, cellosaurus=SlConstants.HEC59_CELLOSAURUS, pmid=pmid)
            sli_list.append(sli)
        msh6 = 'MSH6'
        msh6_sl_partners = {'CKMT2', 'PCK2'}
        for geneB in msh6_sl_partners:
            sli = self.create_sli(geneA=msh6, geneB=geneB,
                                  cell=SlConstants.DLD1_CELL,
                                  cellosaurus=SlConstants.DLD1_CELLOSAURUS, pmid=pmid)
            sli_list.append(sli)
        # HTRA2 NOT synthetic lethal with msh2 or mlh1
        htra2 = 'HTRA2'
        sli = self.create_sli(geneA=msh2, geneB=htra2,
                              cell=SlConstants.HEC59_CELL, cellosaurus=SlConstants.HEC59_CELLOSAURUS,
                              pmid=pmid, SL=False)
        sli_list.append(sli)
        sli = self.create_sli(geneA=mlh1, geneB=htra2,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                              pmid=pmid, SL=False)
        sli_list.append(sli)
        return sli_list

    def martin_2010(self):
        sli_list = []
        mlh1 = 'MLH1'
        polg = 'POLG'
        pmid = "20227038"
        sli = self.create_sli(geneA=mlh1, geneB=polg,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS, pmid=pmid)
        sli_list.append(sli)
        msh2 = 'MSH2'
        polb = 'POLB'
        sli = self.create_sli(geneA=msh2, geneB=polb,
                              cell=SlConstants.HEC59_CELL, cellosaurus=SlConstants.HEC59_CELLOSAURUS, pmid=pmid)
        sli_list.append(sli)
        return sli_list

    def parse(self):
        sli_list = []
        slis = self.martin_2010()
        sli_list.extend(slis)
        slis = self.martin_2011()
        sli_list.extend(slis)
        return sli_list
