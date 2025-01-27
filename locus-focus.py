import argparse
import json
import re
import os
import sys
import time
from zipfile import ZipFile

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez

def readfasta(fp):
	name = None
	seqs = []
	while True:
		line = fp.readline().decode()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

def zreadfasta(zfp, filename):
	with zfp.open(filename) as fp:
		for defline, seq in readfasta(fp):
			m = re.search(r'^(\S+)', defline)
			uid = m.group(1)
			m = re.search(r'organism=([\w ]+)', defline)
			taxname = m.group(1)
			yield taxname, uid, seq

#######
# CLI #
#######

parser = argparse.ArgumentParser()
parser.add_argument('locus', help='json file' )
parser.add_argument('outdir', help='output directory')
parser.add_argument('--clade', default='all',
	help='specifcy clade [%(default)s]')
parser.add_argument('--padding', type=int, default=1000,
	help='exta sequence upstream and downstream [%(default)i]')
parser.add_argument('--delay', type=float, default=0.5,
	help='delay between requests to prevent denial of service')
arg = parser.parse_args()

############
# Download #
############

zfile = f'{arg.locus}-{arg.clade}.zip'
if not os.path.exists(zfile):
	os.system(f'datasets download gene symbol {arg.locus} --ortholog {arg.clade} --filename {zfile}')

###########
# Process #
###########

data = {}
with ZipFile(zfile) as zfp:

	# genomic locations
	with zfp.open('ncbi_dataset/data/data_report.jsonl') as fp:
		for record in fp:
			d = json.loads(record)
			#print(json.dumps(d, indent=2))
			if 'annotations' not in d: continue

			taxid = d['taxId']
			taxname = d['taxname']
			common = d['commonName'] if 'commonName' in d else 'N/A'

			for ann in d['annotations']:
				locs = ann['genomicLocations']
				if len(locs) != 1: sys.exit('unexpected multiple loci')
				loc = locs[0]
				ass = ann['assemblyAccession']
				beg = int(loc['genomicRange']['begin'])
				end = int(loc['genomicRange']['end'])
				ori = loc['genomicRange']['orientation']
				acc = loc['genomicAccessionVersion']
				chm = loc['sequenceName'] if 'sequenceName' in loc else None

				# save
				if taxname not in data:
					data[taxname] = {
						'taxid': taxid,
						'common': common,
						'assembly': ass,
						'acc': acc,
						'chrom': chm,
						'loci': set(),
						'proteins': set(),
						'rnas': set(),
					}
				data[taxname]['loci'].add((acc, beg, end, ori))

	# proteins
	for taxname, uid, seq in zreadfasta(zfp, 'ncbi_dataset/data/protein.faa'):
		data[taxname]['proteins'].add( (uid, seq) )

	# rnas
	for taxname, uid, seq in zreadfasta(zfp, 'ncbi_dataset/data/rna.fna'):
		data[taxname]['rnas'].add( (uid, seq) )

# convert sets to dictionaries
for taxname, record in data.items():
	loci = [{'acc':acc, 'beg':beg, 'end':end, 'ori':ori} for acc, beg, end, ori in data[taxname]['loci']]
	prots = {uid:seq for uid, seq in record['proteins']}
	rnas = {uid:seq for uid, seq in record['rnas']}
	data[taxname]['loci'] = loci
	data[taxname]['proteins'] = prots
	data[taxname]['rnas'] = rnas

# download genome packages for each locus
os.system(f'mkdir -p {arg.outdir}')
for taxname, record in data.items():
	pkg = f'{arg.outdir}/{taxid}.zip'
	if os.path.exists(pkg): continue
	print(json.dumps(record, indent=2))
	taxid = record['taxid']
	chm = record['chrom']
	beg = record['loci'][0]['beg'] - arg.padding
	end = record['loci'][0]['end'] + arg.padding
	location = f'--chromosomes {chm}:{beg}-{end}'
	print(f'downloding {taxid}', file=sys.stderr, flush=True)
	os.system(f'datasets download genome taxon {taxid} {location} --filename {pkg} --include rna,gff3')
	sys.exit()
	time.sleep(arg.delay)

