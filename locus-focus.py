import argparse
import json
import re
import os
import sys
import time
from zipfile import ZipFile

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez

#######
# CLI #
#######

parser = argparse.ArgumentParser()
parser.add_argument('locus', help='json file' )
parser.add_argument('--clade', default='all',
	help='specifcy clade [%(default)s]')
parser.add_argument('--delay', type=float, default=0.5,
	help='delay between requests to prevent denial of service')
arg = parser.parse_args()

zfile = f'{arg.locus}-{arg.clade}.zip'
if not os.path.exists:
	os.system(f'datasets download gene symbol {arg.locus} --ortholog all')
	os.system(f'mv ncbi_dataset.zip {zfile}')

with ZipFile(zfile) as zfp:
	with zfp.open('ncbi_dataset/data/data_report.jsonl') as fp:
		for record in fp:
			d = json.loads(record)
			if 'annotations' not in d: continue

			taxid = d['taxId']
			taxname = d['taxname']
			common = d['commonName'] if 'commonName' in d else 'N/A'

			for ann in d['annotations']:
				locs = ann['genomicLocations']
				if len(locs) != 1: sys.exit('unexpected multiple loci')
				loc = locs[0]
				ass = ann['assemblyAccession']
				beg = loc['genomicRange']['begin']
				end = loc['genomicRange']['end']
				ori = loc['genomicRange']['orientation']
				acc = loc['genomicAccessionVersion']

				print('\t'.join((taxid, taxname, common, ass, acc, beg, end, ori)))


"""
useful files
ncbi_dataset/data/rna.fna
ncbi_dataset/data/protein.faa
ncbi_dataset/data/data_report.jsonl
"""
