from idg2sl.parsers.sl_constants import SlConstants
from idg2sl import SyntheticLethalInteraction


class ManualEntry:
    """
    In a number of papers, only one or a handful of synthetic lethal interactions are described.
    These are very valuable. It is easiest to enter this information by hand.
    """

    def __init__(self):
        self.entries = []
        self._add_reid_2016()
        self._add_parthak_2015()
        self._add_sultana_2013()
        self._add_hoffmann_2014()
        self._add_scholl_2009()
        self._add_barbie_2009()
        self._add_chen_2020()
        self._add_yamada_2020()
        self._add_lelij_2020()
        self._add_liu_2020()
        self._add_huang_2020()
        self._add_szymanska_2020
        self._add_hu_2020()
        self._add_hu_2020a()
        self._add_villalba_2019()
        self._add_paul_2020()
        self._add_chatterjee2019()
        self._add_parvin_2019()
        self._add_park_2019()
        self._add_burdova_2019()
        self._add_zhang2019()
        self._add_li2019()
        self._add_sasaki2019()
        self._add_karpelMassler_2017()
        self._add_li_2018()
        self._add_bajrami_2018()
        self._add_jenkins_2018()
        self._get_gong_2019()

    def create_sli(self, geneA, geneAid, geneB, geneBid, geneApert, geneBpert, cell, cellosaurus, cancer, ncit,
                   assay, pmid,  effecttype = SlConstants.N_A, effectsize=SlConstants.N_A):
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
                                         pmid=pmid,
                                         SL=True)
        return sli


    def _get_gong_2019(self):
        aurka = 'AURKA'
        rb1 = 'RB1'
        pmid = '30373917'
        sli = self.create_sli(geneA=rb1, geneAid=SlConstants.RB1_GENE_ID,geneB=aurka, geneBid=SlConstants.AURKA_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.PHARMACEUTICAL,
                              cell='RB1-Mutant Cancer Cells', cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_jenkins_2018(self):
        tnk2 = 'TNK2'
        ptpn11 = 'PTPN11'
        pmid= '30018082'
        sli = self.create_sli(geneA=ptpn11, geneAid=SlConstants.PTPN11_GENE_ID, geneB=tnk2, geneBid=SlConstants.TNK2_GENE_ID,
                              geneApert=SlConstants.ACTIVATING_MUTATION,geneBpert=SlConstants.PHARMACEUTICAL,
                              cell=SlConstants.N_A,cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.JMML, ncit=SlConstants.JMML_NCIT,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_bajrami_2018(self):
        """
        The cell adhesion glycoprotein E-cadherin (CDH1) is commonly inactivated in breast tumors.
        In order to identify E-cadherin synthetic lethal effects from our MCF7 isogenic cell line siRNA screen,
        we calculated the difference in siRNA Z scores between E-cadherin defective and E-cadherin proficient cells
        and identified 104 E-cadherin synthetic lethal effects (p < 0.05, Fig. 1E; Supplementary Table S3).
        The authors go on to validate ROS1 in detail and do not further analyze these candidates.
        The supplement 3 is an image that cannot easily be ingested -- therefore, we ingest only this highly
        validated SL interaction.
        """
        cdh1 = 'CDH1'
        ros1 = 'ROS1'
        pmid = '29610289'
        sli = self.create_sli(geneA=cdh1, geneAid=SlConstants.CDH1_GENE_ID, geneB=ros1,
                              geneBid=SlConstants.ROS1_GENE_ID,
                              geneApert=SlConstants.CRISPR_CAS9, geneBpert=SlConstants.SI_RNA,
                              cell=SlConstants.MCF7_CELL, cellosaurus=SlConstants.MCF7_CELLOSAURUS,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_li_2018(self):
        atm = 'ATM'
        pten = 'PTEN'
        pmid = '29522753'
        sli = self.create_sli(geneA=pten, geneAid=SlConstants.PTEN_GENE_ID, geneB=atm, geneBid=SlConstants.ATM_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION,geneBpert=SlConstants.PHARMACEUTICAL,
                              cell=SlConstants.MDAMB468_CELL, cellosaurus=SlConstants.MDAMB468_CELLOSAURUS,
                              cancer=SlConstants.N_A,ncit=SlConstants.N_A,
                              assay=SlConstants.CELL_VIABILITY_ASSAY,pmid=pmid)
        self.entries.append(sli)


    def _add_karpelMassler_2017(self):
        """
        bcl-xL is BCL2L1
        """
        bcl2l1 = 'BCL2L1'
        idh1 = 'IDH1'
        pmid = '29057925'
        sli = self.create_sli(geneA=idh1, geneAid=SlConstants.IDH1_GENE_ID,geneB=bcl2l1,geneBid=SlConstants.BCL2L1_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.PHARMACEUTICAL,
                              cell=SlConstants.T98G_CELL, cellosaurus=SlConstants.T98G_CELLOSAURUS,
                              cancer=SlConstants.GLIOBLASTOMA, ncit=SlConstants.GLIOBLASTOMA_NCIT,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY,pmid=pmid)
        self.entries.append(sli)




    def _add_sasaki2019(self):
        parg = 'PARG'
        dusp22 = 'DUSP22'
        pmid = '31142510'
        sli = self.create_sli(geneA=parg, geneAid=SlConstants.PARG_GENE_ID,
                              geneB=dusp22, geneBid=SlConstants.DUSP22_GENE_ID,
                              geneApert=SlConstants.SI_RNA, geneBpert=SlConstants.SH_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.A549_CELL, cellosaurus=SlConstants.A549_CELLOSAURUS,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.SH_RNA,pmid=pmid)
        self.entries.append(sli)

    def _add_li2019(self):
        """
        Note CBP is now called CREBBP
        p300 is now EP300
        """
        crebbp = 'CREBBP'
        ep300 = 'EP300'
        pmid = '31223286'
        sli = self.create_sli(geneA=crebbp, geneAid=SlConstants.CREBBP_GENE_ID,geneB=ep300, geneBid=SlConstants.EP300_GENE_ID,
                              geneApert=SlConstants.SI_RNA, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.CELL_5637, cellosaurus=SlConstants.CELL_5637_CELLOSAURUS,
                              cancer=SlConstants.BLADDER_CARCINOMA, ncit=SlConstants.BLADDER_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_zhang2019(self):
        """
        double inhibition of Aur-A and mTOR showed significant synergistic effects in TNBC cell lines and a xenograft model,
        """
        mtor = 'MTOR'
        aurka = 'AURKA'
        pmid = '31406104'
        sli = self.create_sli(geneA=mtor, geneAid=SlConstants.MTOR_GENE_ID, geneB=aurka, geneBid=SlConstants.AURKA_GENE_ID,
                              geneApert=SlConstants.PHARMACEUTICAL,
                              geneBpert=SlConstants.PHARMACEUTICAL,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell='TNBC cell lines', cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.TRIPLE_NEG_BREAST_CARCINOMA, ncit=SlConstants.TRIPLE_NEG_BREAST_CARCINOMA_NCIT,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY,pmid=pmid)
        self.entries.append(sli)


    def _add_burdova_2019(self):
        ccnf = 'CCNF' # cyclin F
        chek1 = 'CHEK1' # the authors use chk1
        pmid = '31424118'
        sli = self.create_sli(geneA=chek1, geneAid=SlConstants.CHEK1_GENE_ID,
                              geneB=ccnf, geneBid=SlConstants.CCNF_GENE_ID,
                              geneApert=SlConstants.PHARMACEUTICAL,
                              geneBpert=SlConstants.CRISPR_CAS9,
                              effecttype=SlConstants.N_A,
                              effectsize=SlConstants.N_A,
                              cell='cyclin F knock‐out cell line', cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.CELL_VIABILITY_ASSAY,pmid=pmid)
        self.entries.append(sli)

    def _add_park_2019(self):
        myc = 'MCY'
        aurka = 'AURKA'
        pmid = '31429028'
        sli = self.create_sli(geneA=myc, geneAid=SlConstants.MYC_GENE_ID, geneB=aurka,
                              geneBid=SlConstants.AURKA_GENE_ID,
                              geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.PHARMACEUTICAL,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell='Myc-overexpressing lymphoma cell lines', cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.PATIENT_DERIVED_XENOGRAFT, pmid=pmid)
        self.entries.append(sli)

    def _add_parvin_2019(self):
        lmo2 = 'LMO2'
        parp1 = 'PARP1'
        pmid = '31447348'
        sli = self.create_sli(geneA=lmo2, geneAid=SlConstants.LMO2_GENE_ID, geneB=parp1,
                              geneBid=SlConstants.PARP1_GENE_ID,
                              geneApert=SlConstants.OVEREXPRESSION, geneBpert=SlConstants.PHARMACEUTICAL,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.DOHH2_CELL, cellosaurus=SlConstants.DOHH2_CELLOSAURUS,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_chatterjee2019(self):
        pten = 'PTEN'
        pdk1 = 'PDK1'  # Note the authors use the wrong symbol ('PDHK1') for pyruvate dehydrogenase kinase 1
        pmid = '31461649'
        sli = self.create_sli(geneA=pten, geneAid=SlConstants.PTEN_GENE_ID, geneB=pdk1,
                              geneBid=SlConstants.PDK1_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.SH_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.N_A, cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_reid_2016(self):
        """
         PLK1 NCBI Gene id 5347 and CKS1B NCBI Gene id 1163 
         """
        pmid = '27558135'
        cell_line = 'multiple.breast.cancer.cell.lines'
        sli = self.create_sli(geneA='PLK1', geneAid=SlConstants.PLK1_GENE_ID,
                              geneB="CKS1B", geneBid=SlConstants.CKS1B_GENE_ID,
                              geneApert=SlConstants.SH_RNA, geneBpert=SlConstants.OVEREXPRESSION,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=cell_line, cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.BREAST_CARCINOMA, ncit=SlConstants.BREAST_CARCINOMA_NCIT,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                              pmid=pmid)
        self.entries.append(sli)

    def _add_parthak_2015(self):
        """
        # just 2 SL Interactions, hardcoded
        # SRC Gene is proto-oncogene blocked by Dasatinib
        # trying to maximise the Dasatinib sensitivity by SL interaction
        """
        src = 'SRC'
        csnk2a1 = "CSNK2A1"
        CSNK2A1_GENE_ID = "NCBIGene:1457"
        pmid = '26437225'
        effect_type = "correlation"
        sli = self.create_sli(geneA=src, geneAid=SlConstants.SRC_GENE_ID,
                              geneB=csnk2a1, geneBid=CSNK2A1_GENE_ID,
                              geneApert=SlConstants.PHARMACEUTICAL, geneBpert=SlConstants.COHORT_STUDY,
                              effecttype=effect_type, effectsize='-0.82',
                              cell=SlConstants.N_A, cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.RECURRENT_OVARIAN_CANCER,
                              ncit=SlConstants.RECURRENT_OVARIAN_CANCER_NCIT,
                              assay=SlConstants.PHARAMACEUTICAL_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)
        prkce = "PRKCE"
        PRKCE_GENE_ID = "NCBIGene:5581"
        sli = self.create_sli(geneA=src, geneAid=SlConstants.SRC_GENE_ID,
                              geneB=prkce, geneBid=PRKCE_GENE_ID,
                              geneApert=SlConstants.PHARMACEUTICAL, geneBpert=SlConstants.COHORT_STUDY,
                              effecttype=effect_type, effectsize='-0.96',
                              cell=SlConstants.N_A, cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.RECURRENT_OVARIAN_CANCER,
                              ncit=SlConstants.RECURRENT_OVARIAN_CANCER_NCIT,
                              assay=SlConstants.PHARAMACEUTICAL_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_sultana_2013(self):
        atr = 'ATR'
        atr_id = 'NCBIGene:545'
        xrcc1 = 'XRCC1'
        xrcc1_id = 'NCBIGene:7515'
        pmid = '23451157'
        sli = self.create_sli(geneA=atr, geneAid=atr_id, geneB=xrcc1, geneBid=xrcc1_id,
                              geneApert=SlConstants.PHARMACEUTICAL, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.OVCAR3_CELL, cellosaurus=SlConstants.OVCAR3_CELLOSAURUS,
                              cancer=SlConstants.OVARIAN_SEROUS_ADENOCARCINOMA,
                              ncit=SlConstants.OVARIAN_SEROUS_ADENOCARCINOMA_NCIT,
                              assay=SlConstants.CISPLATIN_CYTOTOXICITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_hoffmann_2014(self):
        """
        BRG1 and BRM identified as a SL pair in various different cancer cells
        BRM and BRG1 are closely related paralogs that function as mutually exclusive ATP-dependent catalytic subunits
        of the mSWI/SNF complexes, including A-549 (CVCL_0023)
        The current symbol for BRM is SMARCA2
        The current symbol for BRG1 is SMARCA4

        """
        smarca2 = 'SMARCA2'
        smarca2_id = 'NCBIGene:6595'
        smarca2_perturbation = SlConstants.SI_RNA
        smarca4 = 'SMARCA4'
        smarca4_id = 'NCBIGene:6597'
        smarca4_perturbation = SlConstants.LOF_MUTATION
        pmid = '31427792'
        sli = self.create_sli(geneA=smarca2, geneAid=smarca2_id, geneB=smarca4, geneBid=smarca4_id,
                              geneApert=smarca2_perturbation, geneBpert=smarca4_perturbation,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.A549_CELL, cellosaurus=SlConstants.A549_CELLOSAURUS,
                              cancer=SlConstants.LUNG_ADENOCARCINOMA, ncit=SlConstants.LUNG_ADENOCARCINOMA_NCIT,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_barbie_2009(self):
        """
        Paper is about one SLI
        """
        kras = 'KRAS'
        tbk1 = 'TBK1'
        tbk1_id = 'NCBIGene:29110'
        pmid = '19847166'
        sli = self.create_sli(geneA=kras, geneAid=SlConstants.KRAS_GENE_ID, geneB=tbk1, geneBid=tbk1_id,
                              geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.N_A, cellosaurus=SlConstants.N_A,
                              assay=SlConstants.RNA_INTERFERENCE_ASSAY,
                              cancer=SlConstants.LUNG_ADENOCARCINOMA, ncit=SlConstants.LUNG_ADENOCARCINOMA_NCIT,
                              pmid=pmid)
        self.entries.append(sli)

    def _add_scholl_2009(self):
        """
        KRAS and  STK33
        Scholl C, Fröhling S, Dunn IF, et al. Synthetic lethal interaction between oncogenic KRAS dependency and
        STK33 suppression in human cancer cells. Cell. 2009;137(5):821-834.
        """
        kras = 'KRAS'
        kras_id = 'NCBIGene:3845'
        stk33 = 'STK33'
        stk33_id = 'NCBIGene:65975'
        pmid = '19490892'
        sli = self.create_sli(geneA=kras, geneAid=SlConstants.KRAS_GENE_ID, geneB=stk33, geneBid=stk33_id,
                              geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_chen_2020(self):
        pten = 'PTEN'
        mcl1 = 'MCL1'
        mcl1_id = 'NCBIGene:4170'
        pmid = '32737157'
        cell_line = 'isogeneic GBM cell lines'
        sli = self.create_sli(geneA=pten, geneAid=SlConstants.PTEN_GENE_ID, geneB=mcl1, geneBid=mcl1_id,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.PHARMACEUTICAL,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=cell_line, cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.PHARAMACEUTICAL_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_yamada_2020(self):
        arid1a = 'ARID1A'
        arid1a_id = SlConstants.ARID1A_GENE_ID
        ezh2 = 'EZH2'
        ezh2id = SlConstants.EZH2_GENE_ID
        pmid = '32506298'
        cell_line = 'isogeneic GBM cell lines'
        sli = self.create_sli(geneA=arid1a, geneAid=arid1a_id, geneB=ezh2, geneBid=ezh2id,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.PHARMACEUTICAL,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=cell_line, cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.GASTRIC_CARCINOMA, ncit=SlConstants.GASTRIC_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_lelij_2020(self):
        stag1 = 'STAG1'
        stag1id = 'NCBIGene:10274'
        stag2 = 'STAG2'
        stag2id = 'NCBIGene:10735'
        pmid = '32467316'
        sli = self.create_sli(geneA=stag1, geneAid=stag1id, geneB=stag2, geneBid=stag2id,
                              geneApert=SlConstants.CRISPR_CAS9, geneBpert=SlConstants.LOF_MUTATION,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.KBM7_CELL, cellosaurus=SlConstants.KBM7_CELLOSAURUS,
                              cancer=SlConstants.CML_BCRABL_POS, ncit=SlConstants.CML_BCRABL_POS_NCIT,
                              assay=SlConstants.CRISPR_CAS9_INTERFERENCE_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_liu_2020(self):
        """
        Silencing p300 or CBP selectively inhibited the viability of PTEN-/- HCT116 cells, but not PTEN+/+ ones,
        indicating that the inhibition of p300 and CBP was likely to mediate the synthetic lethality phenotype induced
        by AA
        """
        pten = 'PTEN'
        ep300 = 'EP300'
        pmid = "32398948"
        sli = self.create_sli(geneA=pten, geneAid=SlConstants.PTEN_GENE_ID,
                              geneB=ep300, geneBid=SlConstants.EP300_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                              cancer=SlConstants.COLON_CARCINOMA, ncit=SlConstants.COLON_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)
        crebbp = 'CREBBP'
        sli = self.create_sli(geneA=crebbp, geneAid=SlConstants.CREBBP_GENE_ID,
                              geneB=ep300, geneBid=SlConstants.EP300_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                              cancer=SlConstants.COLON_CARCINOMA, ncit=SlConstants.COLON_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_huang_2020(self):
        ezh2 = 'EZH2'
        pbrm1 = 'PBRM1'
        pmid = '32093567'
        sli = self.create_sli(geneA=ezh2, geneAid=SlConstants.EZH2_GENE_ID,
                              geneB=pbrm1, geneBid=SlConstants.PBRM1_GENE_ID,
                              geneApert=SlConstants.PHARMACEUTICAL, geneBpert=SlConstants.LOF_MUTATION,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.A704_CELL, cellosaurus=SlConstants.A704_CELLOSAURUS,
                              cancer=SlConstants.RENAL_CELL_CARCINOMA, ncit=SlConstants.RENAL_CELL_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_szymanska_2020(self):
        vps4a = 'VPS4A'
        vps4a_id = 'NCBIGene:27183'
        vps4b = 'VPS4B'
        vps4b_id = 'NCBIGene:9525'
        pmid = '31930723'
        sli = self.create_sli(geneA=vps4a, geneAid=vps4a_id, geneB=vps4b, geneBid=vps4b_id,
                              geneApert=SlConstants.SI_RNA, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                              cancer=SlConstants.COLON_CARCINOMA, ncit=SlConstants.COLON_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY,
                              pmid=pmid)
        self.entries.append(sli)

    def _add_hu_2020(self):
        kras = 'KRAS'
        SLC7A11 = 'SLC7A11'
        pmid = '31874110'
        sli = self.create_sli(geneA=kras, geneAid=SlConstants.KRAS_GENE_ID,
                              geneB=SLC7A11, geneBid=SlConstants.SLC7A11_GENE_ID,
                              geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell='KRAS isogenic cell lines', cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.LUNG_ADENOCARCINOMA, ncit=SlConstants.LUNG_ADENOCARCINOMA_NCIT,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_hu_2020a(self):
        tbk1 = 'TBK1'
        vhl = 'VHL'
        pmid = '31810986'
        sli = self.create_sli(geneA=vhl, geneAid=SlConstants.VHL_GENE_ID, geneB=tbk1, geneBid=SlConstants.TBK1_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.UMRC6_CELL, cellosaurus=SlConstants.UMRC6_CELLOSAURUS,
                              cancer=SlConstants.RENAL_CELL_CARCINOMA, ncit=SlConstants.RENAL_CELL_CARCINOMA_NCIT,
                              assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_villalba_2019(self):
        tmprss4 = 'TMPRSS4'
        ddr1 = 'DDR1'
        pmid = '31659178'
        sli = self.create_sli(geneA=tmprss4, geneAid=SlConstants.TMPRSS4_GENE_ID,
                              geneB=ddr1, geneBid=SlConstants.DDR1_GENE_ID,
                              geneApert=SlConstants.SI_RNA,
                              geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.H358_CELL, cellosaurus=SlConstants.H358_CELLOSAURUS,
                              cancer=SlConstants.MINIMALLY_INVASIVE_LUNG_ADENOCARCINOMA,
                              ncit=SlConstants.MINIMALLY_INVASIVE_LUNG_ADENOCARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_paul_2020(self):
        """
        G - Quadruplex - Binding  Small Molecule Induces  Synthetic  Lethality in Breast
        Cancer Cells by Inhibiting c - MYC and BCL2 Expression
        """
        myc = "MYC"
        bcl2 = "BCL2"

        pmid = "31621996"
        sli = self.create_sli(geneA=myc, geneAid=SlConstants.MYC_GENE_ID, geneB=bcl2, geneBid=SlConstants.BCL2_GENE_ID,
                              geneApert=SlConstants.SI_RNA, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.MCF7_CELL, cellosaurus=SlConstants.MCF7_CELLOSAURUS,
                              cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                              assay=SlConstants.CYTOTOXICITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def get_entries(self):
        return self.entries
