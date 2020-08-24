from idg2sl import *
import os

## First download (if needed) and parse the HGNC file with symbol/NCBI Gene/Ensembl mappings
hgnc_fname = SL_DatasetParser.get_local_hgncfile_name()
if not os.path.exists(hgnc_fname):
    SL_DatasetParser.get_hgnc_file()

hgnc = HgncParser(hgnc_fname)
entrez_dict = hgnc.get_entrez_dictionary()
ensembl_dict = hgnc.get_ensembl_dictionary()
synonym_dict = hgnc.get_synonym_dictionary()


def show_stats(name, sli_list):
    pos = sum(sl.get_SL() for sl in sli_list)
    neg = sum(not sl.get_SL() for sl in sli_list)
    print("[INFO] %s: %d positive and %d negative entries" % (name, pos, neg))


baldwin2010 = Baldwin2010Parser()
baldwin2010_list = baldwin2010.parse()
show_stats("Baldwin et al 2010", baldwin2010_list)


blomen2015 = Blomen2015Parser()
blomen2015_list = blomen2015.parse()
show_stats("Blomen et al 2015", blomen2015_list)

bommi2008 = Bommi2008Parser()
bommi2008_list = bommi2008.parse()
show_stats("Bommi et al 2008", bommi2008_list)

brough2018 = Brough2018Parser()
brough2018_list = brough2018.parse()
show_stats("Brough et al 2018", brough2018_list)

chakraborty2017 = Chakraborty2017Parser()
chakraborty2017_list = chakraborty2017.parse()
show_stats("chakraborty et al 2017", chakraborty2017_list)

dai2013 = Dai2013Parser()
dai2013_list = dai2013.parse()
show_stats("Dai et al 2013", dai2013_list)

han2017 = Han2017Parser()
han2017_list = han2017.parse()
show_stats("Han et al 2017", han2017_list)

josse2014 = Josse2014Parser()
josse2014_list = josse2014.parse()
show_stats("Josse et al 2014", josse2014_list)

kang2015 = Kang2015Parser()
kang2015_list = kang2015.parse()
show_stats("Kang et al 2015", kang2015_list)

kessler2012 = Kessler2012Parser()
kessler2012_list = kessler2012.parse()
show_stats("Kessler et al 2012", kessler2012_list)

krastev2011 = Krastev2011Parser()
krastev2011_list = krastev2011.parse()
show_stats("Krastev et al 2011", krastev2011_list)

lord2008 = Lord2008Parser()
lord2008_list = lord2008.parse()
show_stats("Lord et al 2008", lord2008_list)

# Luo et al 2009
luo2009parser = Luo2009Parser()
luo2009_list = luo2009parser.parse()
show_stats("Luo et al 2009", luo2009_list)

manual = ManualEntry(entrez=entrez_dict, ensembl=ensembl_dict, synonym=synonym_dict)
manual_list = manual.get_entries()
show_stats("Manually entered single-SLI studies (part zero)", manual_list)

manual_one = ManualEntryOne(entrez=entrez_dict, ensembl=ensembl_dict, synonym=synonym_dict)
manual_one_list = manual_one.get_entries()
show_stats("Manually entered single-SLI studies (part one)", manual_one_list)

martin2010 = Martin2010and2011Parser()
martin2010_list = martin2010.parse()
show_stats("Martin et al 2010/2011", martin2010_list)


mengwasser_2019 = Mengwasser2019Parser()
mengwasser_2019_list = mengwasser_2019.parse()
show_stats("Mengwasser et al 2019", mengwasser_2019_list)

mohni2014 = Mohni2014Parser(entrez=entrez_dict, ensembl=ensembl_dict, synonym=synonym_dict)
mohni2014_list = mohni2014.parse()
show_stats("Mohni et al 2014", mohni2014_list)

mondal2019 = Mondal2019Parser()
mondal2019_list = mondal2019.parse()
show_stats("Mondal et al 2019", mondal2019_list)

oser2019 = Oser2019Parser()
oser2019_list = oser2019.parse()
show_stats("Oser et al 2019", oser2019_list)

schick2019 = Schick2019Parser()
schick2019_list = schick2019.parse()
print("[INFO] Schick et al 2019  n=%d SL interactions" % len(schick2019_list))
show_stats("Schick et al 2019", schick2019_list)

shen2015 = Shen2015Parser()
shen2015_list = shen2015.parse()
show_stats("Shen et al 2015 ", shen2015_list)

shen2017 = Shen2017Parser()
shen2017_list = shen2017.parse()
show_stats("Shen et al 2017", shen2017_list)

srivas2016 = Srivas2016Parser()
srivas2016_list = srivas2016.parse()
show_stats("Srivas et al 2016", srivas2016_list)

steckel2012 = Steckel2012Parser()
steckel2012_list = steckel2012.parse()
show_stats("Steckel et al 2012", steckel2012_list)

sullivan2012 = Sullivan2012Parser()
sullivan2012_list = sullivan2012.parse()
show_stats("Sullivan et al 2012", sullivan2012_list)

sun2019 = Sun2019Parser()
sun2019_list = sun2019.parse()
show_stats("Sun et al 2019", sun2019_list)

toyoshima2008 = Toyoshima2008Parser()
toyoshima2008_list = toyoshima2008.parse()
show_stats("Toyoshima et al 2008", toyoshima2008_list)

turner2008 = Turner2008Parser()
turner2008_list = turner2008.parse()
show_stats("Turner et al 2008", turner2008_list)

vizeacoumar2013 = Vizeacoumar2013Parser()
vizeacoumar2013_list = vizeacoumar2013.parse()
show_stats("Vizeacoumar et al 2013", vizeacoumar2013_list)

wang2016 = Wang2016Parser()
wang2016_list = wang2016.parse()
show_stats("Wang et al 2016", wang2016_list)

wang2017 = Wang2017Parser()
wang2017_list = wang2017.parse()
show_stats("Wang et al 2017", wang2017_list)

wang_2019 = Wang2019Parser()
wang_2019_list = wang_2019.parse()
show_stats("Wang et al 2019", wang_2019_list)

sli_lists = [baldwin2010_list, bommi2008_list, blomen2015_list, brough2018_list, chakraborty2017_list, dai2013_list,
             han2017_list, josse2014_list, kang2015_list, kessler2012_list,
             krastev2011_list, lord2008_list, luo2009_list, martin2010_list, mengwasser_2019_list, mohni2014_list,
             mondal2019_list,
             oser2019_list,
             shen2015_list, shen2017_list, schick2019_list, srivas2016_list, steckel2012_list, sullivan2012_list,
             sun2019_list, toyoshima2008_list, turner2008_list, vizeacoumar2013_list, wang2016_list, wang2017_list,
             wang_2019_list, manual_list, manual_one_list]
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

ensembl_file = "SL_data.tsv"
fh = open(ensembl_file, 'wt')
fh.write(SyntheticLethalInteraction.get_tsv_with_ensembl_header() + "\n")
for sli in all_sli_list:
    try:
        fh.write(sli.get_tsv_line_with_ensembl(ensembl_dict) + "\n")
    except:
        print("Bad for %s (pmid:%s) " % (sli.get_gene_A_symbol(), sli.get_pmid()))
fh.close()
