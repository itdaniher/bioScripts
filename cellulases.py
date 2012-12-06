#! /usr/bin/python
# Ian Daniher and Gregory Edelston
# 2012.12.06
# basic automation for resolving gene IDs to strains

# use the beautifulsoup html parsing library
from bs4 import BeautifulSoup
# use the 'requests' http request library
import requests

# use biopy for its BLAST wrapper & parser
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# target gene names, copypasta'd from excel
cellulases = ['M940CN_100140017', 'M940CN_100140015', 'M940CN_100592510', 'M940CN_100762912', 'M940CN_10007303', 'M940CN_10162736', 'M940CN_10009096', 'M940CN_10305052', 'M940CN_10081933', 'M940CN_10033693', 'M940CN_10033692', 'M940CN_10077684', 'M940CN_10199672', 'M940CN_10082643', 'M940CN_10136363', 'M940CN_10133142', 'M940CN_10278721', 'M940CN_10050671', 'M940CN_10331251', 'M940CN_10039541', 'M940CN_10028752', 'M940CN_10219582', 'M940CN_10069291', 'M940CN_10062201', 'M940CN_10329602', 'M940CN_10778581', 'M940CN_10041532', 'M940CN_10188451', 'M940CN_10331921', 'M940CN_10424941', 'M940CN_10638351', 'M940CN_10757501', 'M940CN_10899401', 'M940CN_10690371', 'M940CN_10711241', 'M940CN_10789471', 'M940CN_10579281', 'M940CN_10536621', 'M940CN_10408631', 'M940CN_10999221', 'M940CN_10843621', 'M940CN_10497351', 'M940CN_10672031', 'M940CN_11010611', 'M940CN_10504071', 'M940CN_10850661', 'M940CN_10891431', 'M940CN_10748931', 'M940CN_10872391', 'M940CN_10511111', 'M940CN_10932341', 'M940CN_10858061', 'M940CN_10837811', 'M940CN_10728041', 'M940CN_10872361', 'M940CN_10867131', 'M940CN_10727351', 'M940CN_10286661', 'M940CN_10957521']


def getAAfromGeneID(id):
	# get gene name, return amino acid sequence in FASTA format as string
	JGIRequest = requests.request("GET", "http://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=MetaGeneDetail&page=genePageMainFaa&taxon_oid=3300000504&data_type=assembled&gene_oid="+str(id))
	soup = BeautifulSoup(JGIRequest.text)
	AA = str(''.join(soup.findAll("div", {"id":"content_other"})[0].pre.text.strip()))
	return AA

def getscaffoldIDfromGeneID(id):
	# get gene name, return scaffold ID and size 
	JGIRequest = requests.request("GET", "http://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=MetaGeneDetail&page=metaGeneDetail&data_type=assembled&taxon_oid=3300000504&gene_oid="+str(id))
	soup = BeautifulSoup(JGIRequest.text)
	scaffoldInfo = [ item for item in soup.findAll("th", {"class":"subhead"}) if item.text == "Scaffold Source"][0].next.next.next.a.text.split()
	return scaffoldInfo[0], int(''.join( [letter for letter in scaffoldInfo[1] if letter in map(str, range(10))]))


def blast(record):
	# get amino acid sequence in FASTA format as string, return match names
	print "Blasting"
	# make a BLAST request, get FASTA string, return XML
	result_handle = NCBIWWW.qblast("blastp", "nr", record, format_type="XML")
	# get XML, return python object
	blast_records = NCBIXML.parse(result_handle)
	blast_record = blast_records.next()

	three_best = []
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			three_best.append((alignment,hsp))
			while len(three_best) > 3:
				worst_e = -1
				for i in three_best:
					if i[1].expect > worst_e:
						worst_e = i[1].expect
						to_remove = i
				three_best.remove(to_remove)

	# clean three_best list for returning; syntactic magic
	return map(lambda x: x[0].title.split('|')[4], three_best)

if __name__ == "__main__":
	# iterate through list of cellulase gene IDs
	for cellulase in cellulases:
		aaseq = getAAfromGeneID(cellulase)
		matches = blast(aaseq)
		print cellulase
		print getscaffoldIDfromGeneID(cellulase)
		print matches
		print ''
