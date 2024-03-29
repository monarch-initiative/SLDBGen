from idg2sl.parsers.sl_constants import SlConstants
from idg2sl import SyntheticLethalInteraction
from .manual_entry import ManualEntry


class ManualEntry2(ManualEntry):
    """
       In a number of papers, only one or a handful of synthetic lethal interactions are described.
       These are very valuable. It is easiest to enter this information by hand.
       This class is basically the same as ManualEntry but split up to limit the size of the
       entries
       """

    def __init__(self, entrez, ensembl, synonym):
        super().__init__(entrez=entrez, ensembl=ensembl, synonym=synonym)
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
        self._add_Xu_2019()
        self._add_oike_2013()
        self._add_bajrami_2012()
        self._add_lee_2018()
        self._add_bitler_2015()
        self._add_sinha_2017()
        self._add_morandell_2013()
        self._add_imai_2014()
        self._add_wu_2018()
        self._add_zhou_2014()
        self._add_sajesh_2015()
        self._add_hocke_2016()
        self._add_paul_2016()
        self._add_christodoulou_2017()
        self._add_wang_2010()
        self._add_romero_2014()
        self._add_bian_2014()
        self._add_dietlein_2014()
        self._add_wittig_blaich_2017()
        self._add_molenaar_2009()
        self._add_bajrami_2014()
        self._add_yi_2015()
        self._add_lee_2019()
        self._add_kwok_2016()
        self._add_konda_2017()
        self._add_tiong_2014()
        self._add_bryant_2005()
        self._add_sullivan_reed_2018()
        self._add_dietlein_2015()

    def _add_dietlein_2015(self):
        pmid = '26140595'
        chek1 = 'CHEK1'  # current symbol for Chk1
        mapkapk2 = 'MAPKAPK2'  # current symbol for MK2
        kras = 'KRAS'
        self.create_and_add_sli(geneA=chek1, geneB=mapkapk2, geneApert=SlConstants.PHARMACEUTICAL,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell='KRAS-Driven Cancer Cells',
                                cellosaurus=SlConstants.N_A, background_dependency_gene_symbol=kras,
                                background_dependency_status=SlConstants.ACTIVATING_MUTATION,
                                assay= SlConstants.PATIENT_DERIVED_XENOGRAFT, pmid=pmid)


    def _add_sullivan_reed_2018(self):
        pmid = '29898385'
        rad52 = 'RAD52'
        parp1 = 'PARP1'
        brca1 = 'BRCA1'
        self.create_and_add_sli(geneA=rad52, geneB=parp1, geneApert=SlConstants.PHARMACEUTICAL,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.MDAMB436_CELL,
                                cellosaurus=SlConstants.MDAMB436_CELLOSAURUS, background_dependency_status=SlConstants.LOF_MUTATION,
                                background_dependency_gene_symbol=brca1, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_voss_2017(self):
        pmid = '28138868'
        parp1 = 'PARP1'
        ddx3 = 'DDX3'
        brca1 = 'BRCA1'
        self.create_and_add_sli(geneA=ddx3, geneB=parp1, geneApert=SlConstants.PHARMACEUTICAL,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.MCF7_CELL,
                                cellosaurus=SlConstants.MCF7_CELLOSAURUS, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                background_dependency_status=SlConstants.WILDTYPE,
                                background_dependency_gene_symbol=brca1,
                                pmid=pmid)

    def _add_bryant_2005(self):
        pmid = '15829966'
        brca2 = 'BRCA2'
        parp1 = 'PARP1'
        self.create_and_add_sli(geneA=brca2, geneB=parp1, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.MDAMB231_CELL,
                                cellosaurus=SlConstants.MDAMB231_CELLOSAURUS, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_tiong_2014(self):
        pmid = '24947187'
        tp53 = 'TP53'
        csnk1e = 'CSNK1E'
        self.create_and_add_sli(geneA=tp53, geneB=csnk1e, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SH_RNA, cell=SlConstants.HCT_116,
                                cellosaurus=SlConstants.HCT_116_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)
        ctnnb1 = 'CTNNB1'
        self.create_and_add_sli(geneA=tp53, geneB=ctnnb1, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SH_RNA, cell=SlConstants.HCT_116,
                                cellosaurus=SlConstants.HCT_116_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_konda_2017(self):
        pmid = '28182330'
        hsbp1 = 'HSPB1'  # heat shock protein family B (small) member 1
        met = 'MET'
        self.create_and_add_sli(geneA=met, geneB=hsbp1, geneApert=SlConstants.OVEREXPRESSION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.EBC1_CELL,
                                cellosaurus=SlConstants.EBC1_CELLOSAURUS, assay=SlConstants.APOPTOSIS_ASSAY,
                                pmid=pmid)
        egfr = 'EGFR'
        self.create_and_add_sli(geneA=egfr, geneB=hsbp1, geneApert=SlConstants.OVEREXPRESSION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.DIFI_CELL,
                                cellosaurus=SlConstants.DIFI_CELLOSAURUS, assay=SlConstants.APOPTOSIS_ASSAY,
                                pmid=pmid)
        braf = 'BRAF'
        self.create_and_add_sli(geneA=braf, geneB=hsbp1, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.COLO205_CELL,
                                cellosaurus=SlConstants.COLO205_CELLOSAURUS, assay=SlConstants.APOPTOSIS_ASSAY,
                                pmid=pmid)

    def _add_kwok_2016(self):
        pmid = '26563132'
        atr = 'ATR'
        atm = 'ATM'
        tp53 = 'TP53'
        self.create_and_add_sli(geneA=tp53, geneB=atr, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.MEC1_CELL,
                                cellosaurus=SlConstants.MEC1_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.create_and_add_sli(geneA=atm, geneB=atr, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.PHARMACEUTICAL, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_lee_2019(self):
        pmid = '30266755'
        myc = 'MYC'
        tnfsf10 = 'TNFSF10'  # current symbol for TRAIL
        self.create_and_add_sli(geneA=myc, geneB=tnfsf10, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.SI_RNA, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_yi_2015(self):
        pmid = '25955731'
        egfr = 'EGFR'
        met = 'MET'
        self.create_and_add_sli(geneA=egfr, geneB=met, geneApert=SlConstants.PHARMACEUTICAL,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.MDAMB231_CELL,
                                cellosaurus=SlConstants.MDAMB231_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_bajrami_2014(self):
        pmid = '24240700'
        parp1 = 'PARP1'
        cdk12 = 'CDK12'
        self.create_and_add_sli(geneA=parp1, geneB=cdk12, geneApert=SlConstants.PHARMACEUTICAL,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.OV90_CELL,
                                cellosaurus=SlConstants.OV90_CELLOSAURUS,
                                cancer=SlConstants.OVARIAN_SEROUS_ADENOCARCINOMA,
                                ncit=SlConstants.OVARIAN_SEROUS_ADENOCARCINOMA_NCIT,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_molenaar_2009(self):
        pmid = '19525400'
        mycn = 'MYCN'
        cdk2 = 'CDK2'
        self.create_and_add_sli(geneA=mycn, geneB=cdk2, geneApert=SlConstants.OVEREXPRESSION,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.IMR32_CELL,
                                cellosaurus=SlConstants.IMR32_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_wittig_blaich_2017(self):
        pmid = '28423600'
        braf = 'BRAF'
        dusp6 = 'DUSP6'
        self.create_and_add_sli(geneA=braf, geneB=dusp6, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.SI_RNA, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                cell=SlConstants.N_A,
                                cellosaurus=SlConstants.N_A,
                                cancer=SlConstants.MELANOMA, ncit=SlConstants.MELANOMA_NCIT,
                                pmid=pmid)

    def _add_dietlein_2014(self):
        """
        mutations in PRKDC were strongly associated with KU60648 resistance
        """
        pmid = '24556366'
        msh3 = 'MSH3'
        prkdc = 'PRKDC'
        self.create_and_add_sli(geneA=msh3, geneB=prkdc, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.HCC44_CELL,
                                cellosaurus=SlConstants.HCC44_CELLOSAURUS,
                                assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_bian_2014(self):
        pmid = '24425774'
        ppp2r1a = 'PPP2R1A'
        mad2l1 = 'MAD2L1'  # current symbol for MAD2
        self.create_and_add_sli(geneA=mad2l1, geneB=ppp2r1a, geneApert=SlConstants.OVEREXPRESSION,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.HELA_CELL,
                                cellosaurus=SlConstants.HELA_CELLOSAURUS, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_romero_2014(self):
        """
        Depletion of BRG1 dramatically impaired (by >95%) cell viability in the MAX-deficient cells. To test whether
        this behavior also took place in cancer cells with wild-type MAX, we depleted BRG1 in a panel of six lung
        cancer cell lines with amplification of either MYC, MYCL, or MYCN. Only a moderate decrease (<25%) in cell
        growth was found in some cells, implying that the depletion of BRG1 was preferentially toxic in
         MAX-deficient cells.
        """
        pmid = '24362264'
        max = 'MAX'
        smarca4 = 'SMARCA4'  # previously BRG1
        self.create_and_add_sli(geneA=max, geneB=smarca4, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.LU134A_CELL,
                                cellosaurus=SlConstants.LU134A_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_wang_2010(self):
        pmid = '20562906'
        kras = 'KRAS'
        snai2 = 'SNAI2'  # previous symbol SNAIL2
        self.create_and_add_sli(geneA=kras, geneB=snai2, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.HCT_116,
                                cellosaurus=SlConstants.HCT_116_CELLOSAURUS, assay=SlConstants.SG_RNA_DEPLETION_ASSAY,
                                pmid=pmid)

    def _add_christodoulou_2017(self):
        pmid = '28415695'
        kras = 'KRAS'
        copb2 = 'COPB2'
        self.create_and_add_sli(geneA=kras, geneB=copb2, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.SI_RNA, cell='Homozygous pancreatic KRAS mutated cell lines',
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_paul_2016(self):
        pmid = '27418135'
        ephb6 = 'EPHB6'
        src = 'SRC'
        self.create_and_add_sli(geneA=ephb6, geneB=src, geneApert=SlConstants.PROMOTER_HYPERMETHYLATION,
                                geneBpert=SlConstants.SH_RNA, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_hocke_2016(self):
        pmid = '26755646'
        atr = 'ATR'
        pold1 = 'POLD1'
        self.create_and_add_sli(geneA=atr, geneB=pold1, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.DLD1_CELL,
                                cellosaurus=SlConstants.DLD1_CELLOSAURUS, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_sajesh_2015(self):
        pmid = '26318585'
        sod1 = 'SOD1'
        sl_genes = {'BLM', 'CHEK2'}
        for geneB in sl_genes:
            self.create_and_add_sli(geneA=sod1, geneB=geneB,
                                    geneApert=SlConstants.LOF_MUTATION,
                                    geneBpert=SlConstants.SI_RNA,
                                    cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                    assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_zhou_2014(self):
        pmid = '25495526'
        myc = 'MYC'
        prkdc = 'PRKDC'
        self.create_and_add_sli(geneA=myc, geneB=prkdc, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.SH_RNA, cell=SlConstants.WI38_CELL,
                                cellosaurus=SlConstants.WI38_CELLOSAURUS, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_wu_2018(self):
        pmid = '30097580'
        arid1a = 'ARID1A'
        aurka = 'AURKA'
        self.create_and_add_sli(geneA=arid1a, geneB=aurka, geneApert=SlConstants.CRISPR_CAS9,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.HCT_116,
                                cellosaurus=SlConstants.HCT_116_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_imai_2014(self):
        pmid = '24378760'
        tp53 = 'TP53'  # specific for R175H
        id1 = 'ID1'
        self.create_and_add_sli(geneA=tp53, geneB=id1, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.SKBR3_CELL,
                                cellosaurus=SlConstants.SKBR3_CELLOSAURUS, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_morandell_2013(self):
        pmid = '24239348'
        tp53 = 'TP53'
        MAPKAPK2 = 'MAPKAPK2'  # current symbol for MK2
        self.create_and_add_sli(geneA=tp53, geneB=MAPKAPK2, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SI_RNA, cell=SlConstants.A549_CELL,
                                cellosaurus=SlConstants.A549_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_sinha_2017(self):
        """
        Mainly a computational paper, but the authors validate one SLI
        """
        pmid = '28561042'
        idh1 = 'IDH1'
        acaca = 'ACACA'
        self.create_and_add_sli(geneA=idh1, geneB=acaca, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.THP1_CELL,
                                cellosaurus=SlConstants.THP1_CELLOSAURUS, assay=SlConstants.CELL_VIABILITY_ASSAY,
                                pmid=pmid)

    def _add_bitler_2015(self):
        """
        Highly specific EZH2 inhibitors (such as GSK126)
        """
        pmid = '25686104'
        arid1a = 'ARID1A'
        ezh2 = 'EZH2'
        self.create_and_add_sli(geneA=arid1a, geneB=ezh2, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.PHARMACEUTICAL, cell=SlConstants.RMG1_CELL,
                                cellosaurus=SlConstants.RMG1_CELLOSAURUS, cancer=SlConstants.OVARIAN_CCC,
                                ncit=SlConstants.OVARIAN_CCC_NCIT, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_lee_2018(self):
        pmid = '30101194'
        tlk2 = 'TLK2'
        chk1 = 'CHEK1'
        # we observed DNA damage signaling in cells treated with the checkpoint inhibitor UCN-01 (Fig. 4A and fig. S4A).
        # Co-depletion of TLK2 substantially enhanced this response
        self.create_and_add_sli(geneA=tlk2, geneB=chk1, geneApert=SlConstants.SI_RNA,
                                geneBpert=SlConstants.PHARMACEUTICAL,
                                cell=SlConstants.MDAMB231_CELL, cellosaurus=SlConstants.MDAMB231_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        # Similar results were obtained in U-2-OS cells treated with two different CHK1 inhibitors,
        # in MDA-MB-231 cells treated with a CHK1 inhibitor (AZD7762; fig. S4C), and upon treatment with ATR inhibitors
        # (AZ20 and ETP-46464; fig. S4, D and E) (26–29).
        atr = 'ATR'
        self.create_and_add_sli(geneA=tlk2, geneB=atr, geneApert=SlConstants.SI_RNA,
                                geneBpert=SlConstants.PHARMACEUTICAL,
                                cell=SlConstants.MDAMB231_CELL, cellosaurus=SlConstants.MDAMB231_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        parp1 = 'PARP1'
        # We found that the addition of olaparib strongly decreased the survival of MDA-MB-231 cells depleted for TLK activity
        self.create_and_add_sli(geneA=tlk2, geneB=parp1, geneApert=SlConstants.SI_RNA,
                                geneBpert=SlConstants.PHARMACEUTICAL,
                                cell=SlConstants.MDAMB231_CELL, cellosaurus=SlConstants.MDAMB231_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_bajrami_2012(self):
        pmid = '22933245'
        parp1 = 'PARP1'
        nampt = 'NAMPT'
        self.create_and_add_sli(geneA=parp1, geneB=nampt, geneApert=SlConstants.PHARMACEUTICAL,
                                geneBpert=SlConstants.SI_RNA, cell='panel of TN models', cellosaurus=SlConstants.N_A,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_oike_2013(self):
        pmid = '23872584'
        smarca4 = 'SMARCA4'  # current symbol for BRG1
        smarca2 = 'SMARCA2'  # current symbol for BRM
        self.create_and_add_sli(geneA=smarca2, geneB=smarca4, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SI_RNA, cell='multiple BRG1-deficient cells',
                                cellosaurus=SlConstants.N_A, assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_Xu_2019(self):
        pmid = '30885978'
        hk1 = 'HK1'
        hk2 = 'HK2'
        self.create_and_add_sli(geneA=hk1, geneB=hk2, geneApert=SlConstants.ANTISENSE_OLIGO,
                                geneBpert=SlConstants.LOF_MUTATION, cell=SlConstants.OPM1_CELL,
                                cellosaurus=SlConstants.OPM1_CELLOSAURUS, assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                pmid=pmid)

    def _add_kim_2015(self):
        """
        We chose the  24 genes sensitized at least three of these cell lines to the effects of a Met targeting antibody.
        (Figure 1c).
        """
        pmid = '24662823'
        met = 'MET'
        sli_symbols = {'SERPINA3', 'PARP1', 'AREG', 'ATP1A2', 'BCL2L1', 'AKT2', 'BCR', 'FOS', 'CDKN2C', 'SRF', 'FGFR3',
                       'CALR', 'CRK', 'INSRR', 'CD247', 'CASP1', 'CTTN', 'CHRNA7', 'ITGB3', 'CASP2', 'CD151', 'CCND2',
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
        sod1 = 'SOD1'
        self.create_and_add_sli(geneA=rad54b, geneB=sod1, geneApert=SlConstants.SI_RNA,
                                geneBpert=SlConstants.LOF_MUTATION,
                                cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

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
                                geneBpert=SlConstants.SH_RNA, cell=SlConstants.HCT_116,
                                cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)

    def _add_ward_2017(self):
        """
        Multiple lines of evidence.
        """
        pmid = '28628639'
        fen1 = 'FEN1'
        mre11a = 'MRE11'  # (current symbol for MRE11A)
        self.create_and_add_sli(geneA=fen1, geneB=mre11a,
                                geneApert=SlConstants.PHARMACEUTICAL, geneBpert=SlConstants.LOF_MUTATION,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        atm = 'ATM'
        self.create_and_add_sli(geneA=fen1, geneB=atm, geneApert=SlConstants.PHARMACEUTICAL,
                                geneBpert=SlConstants.SH_RNA,
                                cell=SlConstants.FADU_CELL, cellosaurus=SlConstants.FADU_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_helming_2014(self):
        pmid = '24562383'
        arid1a = 'ARID1A'
        arid1b = 'ARID1B'
        self.create_and_add_sli(geneA=arid1a, geneB=arid1b, geneApert=SlConstants.LOF_MUTATION,
                                geneBpert=SlConstants.SH_RNA, cell=SlConstants.TOV21G_CELL,
                                cellosaurus=SlConstants.TOV21G_CELLOSAURUS,
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
