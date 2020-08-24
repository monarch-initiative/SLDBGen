from idg2sl.parsers.sl_constants import SlConstants
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser



class ManualEntryOne(SL_DatasetParser):
    """
       In a number of papers, only one or a handful of synthetic lethal interactions are described.
       These are very valuable. It is easiest to enter this information by hand.
       This class is basically the same as ManualEntry but split up to limit the size of the
       entries
       """

    def __init__(self, entrez, ensembl, synonym):
        super().__init__(fname=None, pmid=None, entrez=entrez, ensembl=ensembl, synonym=synonym)
        self.entries = []
        self._add_mcmanus_2009()
        self._add_ward_2017()
        self._add_helming_2014()
        self._add_koundinya_2018()
        self._add_turchick_2019()
        self._add_guppy_2017()
        self._add_sayesh_2013()
        self._add_pourdehnad_2013()
        self._add_ali_2018()
        self._add_kim_2015()

    def create_and_add_sli(self, geneA, geneB, geneApert, geneBpert, assay, pmid,
                           cell=SlConstants.N_A, cellosaurus=SlConstants.N_A,
                           cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                           effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                           background_dependency_status=SlConstants.N_A,
                           background_dependency_gene_symbol=SlConstants.N_A,
                           background_dependency_gene_id=SlConstants.N_A,
                           sl=True):
        geneAid = self.get_ncbigene_curie(geneA)
        geneBid = self.get_ncbigene_curie(geneB)
        sli = SyntheticLethalInteraction(gene_A_symbol=geneA,
                                         gene_A_id=geneAid,
                                         gene_B_symbol=geneB,
                                         gene_B_id=geneBid,
                                         gene_A_pert=geneApert,
                                         gene_B_pert=geneBpert,
                                         effect_type=effecttype,
                                         effect_size=effectsize,
                                         cell_line=cell,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         background_dependency_status=background_dependency_status,
                                         background_dependency_gene_symbol=background_dependency_gene_symbol,
                                         background_dependency_gene_id=background_dependency_gene_id,
                                         pmid=pmid,
                                         SL=sl)
        self.entries.append(sli)


    def _add_kim_2015(self):
        """
        We chose the  24 genes sensitized at least three of these cell lines to the effects of a Met targeting antibody.
        (Figure 1c).
        """
        pmid = '24662823'
        met = 'MET'
        metid = self.get_ncbigene_curie(met)
        sli_symbols = {'SERPINA3', 'PARP1', 'AREG', 'ATP1A2', 'BCL2L1', 'AKT2', 'BCR', 'FOS', 'CDKN2C','SRF','FGFR3',
                       'CALR', 'CRK', 'INSRR','CD247', 'CASP1', 'CTTN', 'CHRNA7', 'ITGB3', 'CASP2', 'CD151', 'CCND2',
                       'CTSD', 'EPB41L2'}
        for geneB in sli_symbols:
            self.create_and_add_sli(geneA=met, geneB=geneB, geneApert=SlConstants.INHIBITORY_ANTIBODY,
                                    geneBpert=SlConstants.SI_RNA,
                                    cell=SlConstants.MKN45_CELL, cellosaurus=SlConstants.MKN45_CELLOSAURUS,
                                    assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_ali_2018(self):
        pmid = '30297533'
        parp1 = 'PARP1'
        xrcc1 = 'XRCC1'
        self.create_and_add_sli(geneA=xrcc1, geneB=parp1, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.MDAMB231_CELL,
                                cellosaurus=SlConstants.MDAMB231_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_pourdehnad_2013(self):
        pmid = '23803853'
        # eukaryotic translation initiation factor 4E (eIF4E) binding protein 1 (4EBP1)
        myc = 'MYC'
        eif4ebp1 = 'EIF4EBP1'
        self.create_and_add_sli(geneA=myc, geneB=eif4ebp1, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.LOF_MUTATION, assay=SlConstants.TRANSGENIC_MOUSE_MODEL,
                                pmid=pmid)

    def _add_sayesh_2013(self):
        pmid = '24002644'
        rad54b = 'RAD54B'
        sod1 = 'RAD54B'
        self.create_and_add_sli(geneA=rad54b, geneB=sod1, geneApert=SlConstants.SI_RNA, geneBpert=SlConstants.LOF_MUTATION,
                                cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY,pmid=pmid)

    def _add_guppy_2017(self):
        pmid = '28462496'
        parp1 = 'PARP1'
        rnf20 = 'RNF20'
        self.create_and_add_sli(geneA=parp1, geneB=rnf20, geneApert=SlConstants.SI_RNA, geneBpert=SlConstants.SI_RNA,
                                cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_turchick_2019(self):
        pmid = '30863489'
        pten = 'PTEN'
        rad51 = 'RAD51'
        self.create_and_add_sli(geneA=pten, geneB=rad51, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.INHIBITORY_ANTIBODY, cell=SlConstants.U251MG_CELL,
                                cellosaurus=SlConstants.U251MG_CELLOSAURUS, cancer=SlConstants.MELANOMA,
                                ncit=SlConstants.MELANOMA_NCIT, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_mcmanus_2009(self):
        pmid = '19218431'
        rad54b = 'RAD54B'  # protein kinase C, delta
        fen1 = 'FEN1'
        # These are mouse cells but other experiments were done with human cells
        # that document sufficiently the effect

        self.create_and_add_sli(geneA=rad54b, geneB=fen1, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.SH_RNA, cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)

    def _add_ward_2017(self):
        """
        Multiple lines of evidence.
        """
        pmid = '28628639'
        fen1 = 'FEN1'
        mre11a = 'MRE11' # (current symbol for MRE11A)
        self.create_and_add_sli(geneA=fen1, geneB=mre11a,
                                geneApert=SlConstants.PHARMACEUTICAL, geneBpert=SlConstants.LOF_MUTATION,
                                assay=SlConstants.CELL_VIABILITY_ASSAY,pmid=pmid)
        atm = 'ATM'
        self.create_and_add_sli(geneA=fen1, geneB=atm, geneApert=SlConstants.PHARMACEUTICAL, geneBpert=SlConstants.SH_RNA,
                                cell=SlConstants.FADU_CELL, cellosaurus=SlConstants.FADU_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_helming_2014(self):
        pmid = '24562383'
        arid1a = 'ARID1A'
        arid1b = 'ARID1B'
        self.create_and_add_sli(geneA=arid1a, geneB=arid1b, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SH_RNA,cell=SlConstants.TOV21G_CELL, cellosaurus=SlConstants.TOV21G_CELLOSAURUS,
                                assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)

    def _add_koundinya_2018(self):
        pmid = '29628435'
        kras = 'KRAS'
        dhodh = 'DHODH'
        self.create_and_add_sli(geneA=kras, geneB=dhodh, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.HCT_116,
                                cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)

    def get_entries(self):
        return self.entries
