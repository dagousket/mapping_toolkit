#!/usr/bin/env python2.7
"""
Converts a given BED file into a GFF file with features from a given gff file.
The name of each element should be the list of its features (cf argument id_as_features in gff_to_bed script)
"""

import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Creates GFF file from a given BED file. Note that the features of the GFF are created based on the ID of the BED file.')
parser.add_argument('--bed_file','-f',required=True,help='Name of the BED file to be converted')
parser.add_argument('--source','-s',required=True,help='Name of the source')
parser.add_argument('--mol_type','-m',required=True,help='Name of the molecular type of elements from BED')
parser.add_argument('--is_bed12',required=False,action='store_true',help='Specify this argument if bed file is bed12 formated and contain blocks')
parser.add_argument('--make_gff3', required = False, action = 'store_true', default = False, help = 'Specify if you want to make the output a proper gff3')


def get_Dictbed(file_bed) :
	# put bed info into lis of dictionaries
	with open(file_bed,'rU') as f :
		ListDict = []
		keys = ['chr','start','end','features','strand']
		for l in f :
			values = [l.rstrip('\r\n').split('\t')[i] for i in [0,1,2,3,5]]
			elementDict = dict(zip(keys,values))
			ListDict.append(elementDict)
	f.closed
	return ListDict


def get_Dictbed12(file_bed) :
	# put bed12 info into lis of dictionaries
	with open(file_bed,'rU') as f :
		ListDict = []
		keys = ['chr','start','features','strand','block','size','starting_block']
		for l in f :
			values = [l.rstrip('\r\n').split('\t')[i] for i in [0,1,3,5,9,10,11]]
			elementDict = dict(zip(keys,values))
			ListDict.append(elementDict)
	f.closed
	return ListDict


def get_gff(file_bed,source,mol_type, make_gff3) :
	Dictbed = get_Dictbed(file_bed)
	# write info from dictionary into new gff file
	with open(os.path.splitext(file_bed)[0]+".gff",'w') as f :
		for d in Dictbed :
			if make_gff3 is True:
				f.write('\t'.join([d['chr'],source,mol_type,str(int(d['start'])+1),d['end'],'.',d['strand'],'.','feature_id=' + d['features']+'\n']))
			else:
				f.write('\t'.join([d['chr'],source,mol_type,str(int(d['start'])+1),d['end'],'.',d['strand'],'.',d['features']+'\n']))
	f.closed


def get_gff_from_bed12(file_bed,source,mol_type, make_gff3) :
	Dictbed = get_Dictbed12(file_bed)
	# write info from dictionary into new gff file
	with open(os.path.splitext(file_bed)[0]+".gff",'w') as f :
		for d in Dictbed :
			sizes = map(int,list(d['size'].rstrip(',').split(',')))
			starts = map(int,list(d['starting_block'].rstrip(',').split(',')))
			for i in range(0,int(d['block'])) :
				if make_gff3 is True:
					f.write('\t'.join([d['chr'],source,mol_type,str(starts[i]+1+int(d['start'])),str(starts[i]+1+int(d['start'])+sizes[i]),'.',d['strand'],'.','feature_id=' + d['features']+'\n']))
				else:
					f.write('\t'.join([d['chr'],source,mol_type,str(starts[i]+1+int(d['start'])),str(starts[i]+1+int(d['start'])+sizes[i]),'.',d['strand'],'.',d['features']+'\n']))
	f.closed


if len(sys.argv) == 1 :
	parser.parse_args(['--h'])
else:
	args = parser.parse_args()
	if args.is_bed12 :
		get_gff_from_bed12(args.bed_file,args.source,args.mol_type, args.make_gff3)
	else :
		get_gff(args.bed_file,args.source,args.mol_type, args.make_gff3)
	print("Created GFF file in working directory")

