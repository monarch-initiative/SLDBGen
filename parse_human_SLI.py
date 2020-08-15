from idg2sl.parsers.blomen_2015_parser import Blomen2015Parser
from idg2sl import *
import os


hgnc_fname = SL_DatasetParser.get_local_hgncfile_name()
if not os.path.exists(hgnc_fname):
    SL_DatasetParser.get_hgnc_file()

hgnc = HgncParser(hgnc_fname)
entrez_dict = hgnc.get_entrez_dictionary()
ensembl_dict = hgnc.get_ensembl_dictionary()
synonym_dict = hgnc.get_synonym_dictionary()


vizeacoumar2013 = Vizeacoumar2013Parser()
vizeacoumar2013_list = vizeacoumar2013.parse()
print("[INFO] Vizeacoumar et al 2013  n= %d SL interactions" % len(vizeacoumar2013_list))

brough2018 = Brough2018Parser()
brough2018_list = brough2018.parse()
print("[INFO] Brough et al 2018  n= %d SL interactions" % len(brough2018_list))


kessler2012 = Kessler2012Parser()
kessler2012_list = kessler2012.parse()
print("[INFO] Kessler et al 2012  n= %d SL interactions" % len(kessler2012_list))



manual = ManualEntry()
manual_list = manual.get_entries()
print("We got %d manually entered entries." % len(manual_list))

mohni2014 = Mohni2014Parser(entrez=entrez_dict, ensembl=ensembl_dict, synonym=synonym_dict)
mohni2014_list = mohni2014.parse()
print("[INFO] Mohni et al 2014  n= %d SL interactions" % len(mohni2014_list))


shen2017 = Shen2017Parser()
shen2017_list = shen2017.parse()
print("[INFO] Shen et al 2017  n= %d SL interactions" % len(shen2017_list))


han2017 = Han2017Parser()
han2017_list = han2017.parse()
print("[INFO] Han et al 2017  n= %d SL interactions" % len(han2017_list))


wang2017 = Wang2017Parser()
wang2017_list = wang2017.parse()
print("[INFO] Wang et al 2017  n= %d SL interactions" % len(wang2017_list))


srivas2016 = Srivas2016Parser()
srivas2016_list = srivas2016.parse()
print("[INFO] Srivas et al 2016  n= %d SL interactions" % len(srivas2016_list))


shen2015 = Shen2015Parser()
shen2015_list = shen2015.parse()
print("[INFO] Shen et al 2015  n= %d SL interactions" % len(shen2015_list))


toyoshima2008 = Toyoshima2008Parser()
toyoshima2008_list = toyoshima2008.parse()
print("[INFO] Toyoshima et al 2008  n= %d SL interactions" % len(toyoshima2008_list))


lord2008 = Lord2008Parser()
lord2008_list = lord2008.parse()
print("[INFO] Lord et al 2008  n= %d SL interactions" % len(lord2008_list))


steckel2012 = Steckel2012Parser()
steckel2012_list = steckel2012.parse()
print("[INFO] Steckel et al 2012  n= %d SL interactions" % len(steckel2012_list))

bommi2008 = Bommi2008Parser()
bommi2008_list = bommi2008.parse()
print("[INFO] Bommi et al 2008  n= %d SL interactions" % len(bommi2008_list))

# Turner 2008
turner2008 = Turner2008Parser()
turner_list = turner2008.parse()
print("[INFO] Turner et al 2008  n= %d SL interactions" % len(turner_list))

# Blomen 2015
blomen2015 = Blomen2015Parser()
blomen_list = blomen2015.parse()
print("[INFO] Blomen et al 2015  n= %d SL interactions" % len(blomen_list))

# Luo et al 2009
luo2009parser = Luo2009Parser()
luo2009_list = luo2009parser.parse()
print("[INFO] Luo et al 2009  n= %d SL interactions" % len(luo2009_list))

schick2019 = Schick2019Parser()
schick2019_list = schick2019.parse()
print("[INFO] Schick et al 2019  n= %d SL interactions" % len(schick2019_list))

sli_lists = [bommi2008_list, brough2018_list, han2017_list, kessler2012_list, lord2008_list, luo2009_list,
             mohni2014_list, shen2015_list, shen2017_list, schick2019_list, srivas2016_list, steckel2012_list,
             toyoshima2008_list, turner_list, wang2017_list, manual_list]
all_sli_list = []
for l in sli_lists:
    all_sli_list.extend(l)

n = 0
n_SL = 0
for sli in all_sli_list:
    n += 1
    if sli.get_SL():
        n_SL += 1
print("We got %d interactions including %d synthetic lethal interactions" % (n, n_SL))


def save_SL_data(path, sli_lists):
    with open(path, 'w') as out_f:
        for sli in all_sli_list:
            if sli.get_SL():
                out_f.write(sli.get_gene_A_symbol() + "\t" + sli.get_gene_B_symbol() + "\t")
                out_f.write(str(sli.get_effect_size()) + "\n")

ensembl_file = "SL_data.tsv"
fh = open(ensembl_file, 'wt')
fh.write(SyntheticLethalInteraction.get_tsv_with_ensembl_header() + "\n")
for sli in all_sli_list:
    try:
        fh.write(sli.get_tsv_line_with_ensembl(ensembl_dict) + "\n")
    except:
        print("Bad for ", sli.get_gene_A_symbol())
fh.close()