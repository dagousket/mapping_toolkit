#!/usr/bin/env python2.7
"""Functions to select specific elements from gff file and convert them into BED6 and/or BED12 format"""

import os
import sys
import argparse
import re
parser = argparse.ArgumentParser(description='Creates BED files from a given GFF file with specific filters. In BED12, groups all the elements of a selected molecular type according to their feature type.')
parser.add_argument('--gff_file','-f',required=True,help='Name of the GFF file to be converted')
parser.add_argument('--bed12','-b12',required=False,action='store_true',help='Creates the corresponding BED12 file')
parser.add_argument('--no_bed6','-nb6',required=False,action='store_false',help='Prevents from creating the corresponding BED6 file')
parser.add_argument('--mol_type','-mt',required=False,type=str,default='exon',help='The molecular type (column 3 of the GFF file) selected for the BED files, default is exon')
parser.add_argument('--feature_type','-ft',required=False,type=str,default = 'Parent',help='The feature type (column 9 of the GFF file) selected for the BED files, default is Parent')
parser.add_argument('--id_as_features','-id',required=False,action='store_true',help='Will set the ID of each element as a string containing all its features')
parser.add_argument('--path','-p',required=False,type=str,default='./',help='The location where BED files will be created, default is current working directory')
parser.add_argument('--name','-n',required=False,type=str,default='noinp',help='The name of the BED files, default is the GFF file name')
parser.add_argument('--verbose','-v',required=False,action='store_true',help='Will outpout in stdout the command arguments and the name of each element raising a warning in consistency check')
parser.add_argument('--discard','-d',required=False,action='store_true',help='Will discard the element raising a warning in strand consistency and overlapping check')
parser.add_argument('--skip_exon_number','-s', required = False, action = 'store_true', help = "If set, the program will skip addiing _# for exon number.")

def get_featureDict(featureInfo):
	"""The object "featureInfo" should be the last field of a tab-delimited GFF file.
	And it really must be a GFF3 file, not a GTF or some such nonsense.
	"""
	featureDict = {}
	atoms = [k.lstrip(' ').rstrip(' ') for k in featureInfo.split(';')]
	for k in atoms:
		index, entry = k.split('=')
		featureDict[index.rstrip(' ').lstrip(' ')] = [j.rstrip(' ').lstrip(' ') for j in entry.split(',')]
	return featureDict


def get_Dict(file_gff,mol_type,feature_type) :
	"""Create a dictionary from the GFF file which contains the infomations of all feature_type of a given mol_type (ie exon, CDS).
	Requires get_featureDict
	key = feature name
	value = dictionary for each mol_type in the feature with chr, start, stop and strand keys
	"""
	with open(file_gff,'rU') as f :
		myDict = {}
		for l in f :
			if l.rstrip('\r\n').split('\t')[2] == mol_type :
				#get feature type in key of featDict
				featDict = get_featureDict(l.rstrip('\r\n').split('\t')[-1])
				#get mol_type info
				keys = ['chr','start','stop','strand']
				values = [l.rstrip('\r\n').split('\t')[i] for i in [0,3,4,6]]
				#add the key name if id_as_feature argument is specified
				if args.id_as_features :
					keys.append('name')
					values.append(l.rstrip('\r\n').split('\t')[-1])
				exonDict = dict(zip(keys, values))
				#create dict key=feattype value=moltypeDict1
				for t in featDict[feature_type] :
					if t in myDict :
						#add the exon
						myDict[t].append(exonDict)
					else :
						#create key and add exon
						myDict[t] = [exonDict]
	f.closed
	return myDict


def consistency_check(myDict) :
	"""Check the consistency of the data from gff.
	Return false if chromosome or strand consistency is not respected.
	Return true if respected and order the list of mol_type dictionary according to start postition.
	"""
	alloverlap = []
	allstrand = []
	allchr = 0

	for t in myDict :
		#chromosome consistency
		if len(set([d['chr'] for d in myDict[t]])) != 1 :
			allchr += 1
			print 'WARNING : Chromosome consistency is not satisfied for '+args.feature_type+' '+t
		#strand consistency
		elif len(set([d['strand'] for d in myDict[t]])) != 1 :
			allstrand.append(t)
		#order exon based on start values
		startkey = [int(k['start']) for k in myDict[t]]
		exonsort=[x for (y,x) in sorted(zip(startkey,myDict[t]))]
		myDict[t]=exonsort
		#check overlapping
		initstop = -1
		overlap = False
		for exon in myDict[t] :
			if int(exon['start']) < int(initstop) :
				overlap = True
				#print str(exon['start'])+' < '+str(initstop)+' in '+str(t)
				#Commented line print the position where exon overlapping is detected
			initstop = exon['stop']
		if overlap :
			alloverlap.append(t)
	#Count and print warnings
	if (alloverlap and args.verbose):
		print '\nWarning : overlapping '+args.mol_type+' in '+args.feature_type+' :\n'+','.join(alloverlap)
	if (allstrand and args.verbose):
		print '\nWarning : Strand consistency is not satisfied for '+args.feature_type+' :\n'+t+','.join(allstrand)
	print '\nTotal of :\n\t'+str(len(allstrand))+' strand inconsistencies\n\t'+str(len(alloverlap))+' '+args.feature_type+' with overlapping '+args.mol_type+'\n\t'+str(allchr)+' chromosome inconsistencies'
	#Discard inconsistent element
	if args.discard :
		nbdisc = len(myDict)
		reject = list(set(allstrand) | set(alloverlap))
		for k in reject :
			myDict.pop(k)
		nbdisc -= len(myDict)
		print '\nDiscarded '+str(nbdisc)+' elements with strand inconsistencies and/or overlapping\n'
	return True


def bed6_generator(file_gff,mol_type,feature_type,path,bedname):
	"""Create a simple BED file in the working directory.
	Requires get_Dict.
	"""
	#uses os module
	myDict = get_Dict(file_gff,mol_type,feature_type)
	if consistency_check(myDict) :
		#only allow to pursue script if check script runs correctly
		with open(path.rstrip('/')+'/'+bedname+".bed6",'w') as f6 : #check
			for t in myDict :
				exnum = 0
				for exon in myDict[t] :
					#add a number to exons
					exnum += 1
					if args.id_as_features :
						name = exon['name']
					else :
						if args.skip_exon_number is False:
							name = '_'.join([t,str(exnum)])
						else:
							name = t
					f6.write('\t'.join([exon['chr'],str(int(exon['start'])-1),exon['stop'], name, '.', exon['strand']+'\n']))
		f6.closed


def bed12_generator(file_gff,mol_type,feature_type,path,bedname,check) :
	"""Create a BED12 file in the working directory.
	Requires get_Dict.
	"""
	#uses os module
	myDict = get_Dict(file_gff,mol_type,feature_type)
	if not check :
		if consistency_check(myDict) :
			check = True
	if check :
		#only allow to pursue script if check script runs correctly
		with open(path.rstrip('/')+'/'+bedname+".bed12",'w') as f12 :
			for t in myDict :
				blockSizes = [int(k['stop'])-int(k['start']) for k in myDict[t]]
				blockStart = [int(k['start'])-int(myDict[t][0]['start']) for k in myDict[t]]
				#Handling size 0 exons appearing in gff file
				blockStart = [blockStart[k] for k in range(len(blockSizes)) if blockSizes[k] > 0]
				blockSizes = [blockSizes[k] for k in range(len(blockSizes)) if blockSizes[k] > 0]
				blockCount = str(len(blockSizes))
				#Handling single size 0 exon in a transcript
				if not blockSizes :
					continue
				#Settinf int as string
				blockStart = ','.join(str(e) for e in blockStart)
				blockSizes = ','.join(str(e) for e in blockSizes)
				# check if correct by :
				# blockStart[-1]+myDict[t][0]['start']+blockSizes[-1] == myDict[t][-1]['stop']
				chromStart = str(int(myDict[t][0]['start'])-1)
				chromEnd = myDict[t][-1]['stop']
				if args.id_as_features :
					#recreate the list of feature without the feat_type from argument
					featD = get_featureDict(myDict[t][0]['name'])
					featD['Name'] = [re.sub(':[0-9]$', '', featD['Name'][0])]
					name = ';'.join([str(k)+'='+','.join(v) for (k,v) in featD.items() if k != 'Parent']+["Parent="+str(t)])
				else :
					name = t
				f12.write('\t'.join([myDict[t][0]['chr'],chromStart,chromEnd,name,'0',myDict[t][0]['strand'],chromStart,chromEnd,'0',blockCount,blockSizes,blockStart+'\n']))
		f12.closed




def main(args) :
	"""Creates BED files from GFF file. Options available to select specific molecular type and BED format."""
	if args.name == 'noinp' :
		bedname = os.path.basename(args.gff_file).replace('.gff','')
	else :
		bedname = args.name
	check = False
	if args.no_bed6 :
		print '\nCreating BED6 with all '+args.mol_type+'s from the gff file'
		bed6_generator(args.gff_file,args.mol_type,args.feature_type,args.path,bedname)
		#Avoid redoing consistency check in bed12 if already done for bed6
		check = True
	if args.bed12 :
		print '\nCreating BED12 with all '+args.mol_type+'s from gff file, grouped according to their '+args.feature_type
		bed12_generator(args.gff_file,args.mol_type,args.feature_type,args.path,bedname,check)


if len(sys.argv) == 1 :
	parser.parse_args(['--h'])
else:
	args = parser.parse_args()
	if args.verbose :
		print 'Running gff_to_bed.py with the arguments :\n'
		for x in vars(args) :
			print '--'+x+' '+str(vars(args)[x])
	main(args)





