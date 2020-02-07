

import logging
import os.path
from _collections import defaultdict
from idg2sl import SyntheticLethalInteraction







def parse_luo2009_supplemental_file_S3(path, symbol2entrezID):
    """
    Luo J, et al A genome-wide RNAi screen identifies multiple synthetic lethal
    interactions with the Ras oncogene. Cell. 2009 May 29;137(5):835-48.
    PubMed PMID: 19490893
    Supplemental Table S3 was copied from the original PDF file and stored to
    a file 'luo2009.tsv' in the data directory. The experiment used cell lines
    with activating (oncogenic) mutations in KRAS and used a first screen with
    relative abundance of shRNAs to identify 368 genes using stringent criteria.
    They tested 320 candidates from the frist screen and found 83 shRNAs targeting
    77 genes to preferentially decreased the viability of the KRAS mutant cells
    compared to WT cells, thus indicating SL. They screened 68 of these candidates
    in a second line and could confirm 50 of them (73.5%), indicating that the
    majority of the candidates were likely to be true positives.
    """
    # because of the experiment, geneA is always KRAS. GeneB is in Table S3
    kras_symbol = 'KRAS'
    kras_id = 'NCBIGene:3845'
    kras_perturbation = 'oncogenic_mutation'
    gene2_perturbation = 'shRNA'
    pmid = 'PMID:19490893'
    assays = ['competitive hybridization', 'multicolor competition assay']
    assay_string = ";".join(assays)
    effect_type = 'stddev'

    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Luo et al 2009")
    sli_dict = defaultdict(SyntheticLethalInteraction)

    with open(path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                logging.error("Only got %d fields but was expecting 8" % len(fields))
                i = 0
                for f in fields:
                    print("%d) %s" % (i, f))

                raise ValueError("Malformed line, must have 8 tab-separated fields")

            geneB_sym = fields[0]
            geneB_refSeq = fields[1]
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'

            stddev = float(fields[5])
            SL = True # All data in this set is True # TODO CHECK
            sli = SyntheticLethalInteraction(gene_A_symbol=kras_symbol,
                                             gene_A_id=kras_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=kras_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=stddev,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=SL)
            if geneB_sym in sli_dict:
                # get the entry with the strongest effect size
                sli_b = sli_dict.get(geneB_sym)
                if abs(stddev) > abs(sli_b.effect_size):
                    sli_dict[geneB_sym] = sli
            else:
                # first entry for geneB
                sli_dict[geneB_sym] = sli
    return sli_dict.values()




def parse_bommi_reddi_2008(path, symbol2entrezID):
    """
     Bommi-Reddy A, et al  Kinase requirements in human cells: III. Altered kinase
    requirements in VHL-/- cancer cells detected in a pilot synthetic lethal screen.
    Proc Natl Acad Sci U S A. 2008 Oct 28;105(43):16484-9. PubMed PMID: 18948595.
    The paper describes the use of two cell types: 786-O (Table S4) and RCC4 (Table S5).
    Each of these cells has a loss of function mutation in the tumor suppressor gene VHL.
    The defects in the other gene were induced with a lentiviral vector that produced a shRNA. The
    final results are in Tables S4 and S5.
    The authors require a differential viability above 35% for the 786-O cells or a differential
    viability above 20% for the RCC4 cells. We have copied that values from the supplemental tables
    (which were PDF files) into a TSV file called data/bommi-reddy-2008.tsv. For the output file,
    we demand that at least one shRNA was above threshold, and just report the best outcome per gene.
    The authors showed experimentally that two genes, IRR and HER4, actually are not SL, so we remove these
    by hand
    :return:
    """
    # because of the experiment, geneA is always KRAS. GeneB is in Table S3
    vhl_symbol = 'VHL'
    vhl_id = 'NCBIGene:7428'
    vhl_perturbation = 'lof_mutation'
    gene2_perturbation = 'shRNA'
    pmid = 'PMID:18948595'
    assays = ['competitive hybridization', 'multicolor competition assay']
    assay_string = ";".join(assays)
    effect_type = 'differential_viability'
    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Bommi-Reddy A, et al 2008")
    sli_dict = defaultdict(SyntheticLethalInteraction)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            print(len(fields), line)
            if len(fields) < 4:
                raise ValueError("Only got %d fields but was expecting 4" % len(fields))
            geneB = fields[0]
            if geneB is "IRR" or geneB is "HER4":
                continue
            if geneB in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB))
            else:
                geneB_id = "n/a"
            effect = float(fields[1])
            cell = fields[2]
            table = fields[3]
            assay_string = "differential viability assay {}({})".format(cell, table)
            SL = True  # All data in this set is True # TODO CHECK
            sli = SyntheticLethalInteraction(gene_A_symbol=vhl_symbol,
                                             gene_A_id=vhl_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=vhl_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=effect,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=SL)
            if geneB in sli_dict:
                # get the entry with the strongest effect size
                sli_b = sli_dict.get(geneB)
                if abs(effect) > abs(sli_b.effect_size):
                    sli_dict[geneB] = sli
            else:
                # first entry for geneB
                sli_dict[geneB] = sli
    return sli_dict.values()




def parse_turner_2008(path, symbol2entrezID):
    """
    Parse data from  Turner NC, et al., A synthetic lethal siRNA screen identifying genes mediating
    sensitivity to a PARP inhibitor. EMBO J. 2008 May 7;27(9):1368-77. PubMed PMID: 18388863
    PARP1 was inhibited by the PARP inhibitor KU0058948, and a short interfering RNA library targeting
    779 human protein kinase and kinase assocaited genes was applied.
    :param path:
    :return:
    """
    parp1_symbol = 'PARP1'
    parp1_id = 'NCBIGene:142'
    parp1_perturbation = 'drug'
    gene2_perturbation = 'siRNA'
    pmid = 'PMID:18388863'
    assays = ['competitive hybridization', 'multicolor competition assay']
    assay_string = ";".join(assays)
    effect_type = 'stddev'
    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Turner et al 2008")
    sli_dict = defaultdict(SyntheticLethalInteraction)
    with open(path) as f:
        next(f) # skip header
        for line in f:
            if len(line) < 3:
                raise ValueError("Bad line for Turner et al")
            fields = line.rstrip('\n').split('\t')
            geneB = fields[0]
            if geneB in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB))
            else:
                geneB_id = "n/a"
            zscore = float(fields[1])
            percent_inhib = fields[2]
            if zscore <= -3.0:
                SL = True
            else:
                SL = False
            sli = SyntheticLethalInteraction(gene_A_symbol=parp1_symbol,
                                             gene_A_id=parp1_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=parp1_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=zscore,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=SL)
            if geneB in sli_dict:
                # get the entry with the strongest effect size
                sli_b = sli_dict.get(geneB)
                if abs(zscore) > abs(sli_b.effect_size):
                    sli_dict[geneB] = sli
            else:
                # first entry for geneB
                sli_dict[geneB] = sli
        return sli_dict.values()
