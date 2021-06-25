#!/usr/bin/env python3

import argparse
import copy



def quant2list(quantFile, algorithm):
	accession_Readcounts = list()
	if algorithm == 'agamemnon':
		column = 2
	elif algorithm == 'kallisto':
		column = 3
	with open(quantFile) as f:
		next(f)
		for line in f:
			accession_Readcounts.append([line.split('\t')[0].strip(),\
								float(line.split('\t')[column].strip())])

	return accession_Readcounts



def tempaccessions2taxid(accession_Readcounts, atFile):
	taxid_Accessions = dict()
	not_Found = list()
	accessions = [element[0] for element in accession_Readcounts]
	with open(atFile) as f:
		for line in f:
			if line.split('\t')[1].strip() in taxid_Accessions:
				taxid_Accessions[line.split('\t')[1].strip()].append(line.split('\t')[0].strip())
			else:
				taxid_Accessions[line.split('\t')[1].strip()] = [line.split('\t')[0].strip()]
	
	accessions_Found = set([item for sublist in [key for key in taxid_Accessions.values()] for item in sublist])
	notFound = list(set(accessions) - accessions_Found)

	return (taxid_Accessions, notFound)



def fixEmptyranks(lineage):
	empty = 'NA'
	ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	for key in lineage:
		ranks_Found = list()
		if len(lineage[key]) == 8:
			continue
		else:
			for element in range(0, len(lineage[key])):
				ranks_Found.append(lineage[key][element].split('|')[0])
			for i in range(0, len(ranks)):
				if ranks[i] not in ranks_Found:
					if ranks[i] == 'species':
						lineage[key].insert(1, empty)
					elif ranks[i] == 'genus':
						lineage[key].insert(2, empty)
					elif ranks[i] == 'family':
						lineage[key].insert(3, empty)
					elif ranks[i] == 'order':
						lineage[key].insert(4, empty)
					elif ranks[i] == 'class':
						lineage[key].insert(5, empty)
					elif ranks[i] == 'phylum':
						lineage[key].insert(6, empty)
					elif ranks[i] == 'superkingdom':
						lineage[key].insert(7, empty)

	return lineage



def buildTaxonomy(nodesFile, namesFile, taxid_Accessions):
	nodes = dict()
	names = dict()
	lineage = dict()
	ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
	with open(nodesFile) as f:
		for line in f:
			elements = line.split('|')
			nodes[elements[0].strip()] = [elements[1].strip(), elements[2].strip()]
	with open(namesFile) as f:
		for line in f:
			elements = line.split('|')
			if elements[3].strip() == 'scientific name':
				names[elements[0].strip()] = elements[1].strip()
	for key in taxid_Accessions:
		tempLineage = list()
		helpful_list = list()
		parent = key
		cnt = 0
		while parent not in ['131567', '10239']:
			tempLineage.append(nodes[parent][1] + '|' + names[parent] + '|' + parent)
			parent = nodes[parent][0]
		rank_counter = 0
		for rank in range(0, len(tempLineage)):
			if rank_counter == 0:
				helpful_list.append(tempLineage[rank])
			elif rank_counter > 0:
				if tempLineage[rank].split('|')[0] in ranks:
					helpful_list.append(tempLineage[rank])
			rank_counter += 1
		lineage[key] = helpful_list
	copy_dict = copy.deepcopy(lineage)
	lineage = copy.deepcopy(copy_dict)

	lineage = fixEmptyranks(lineage)

	return lineage



def counts2lineage(accession_Readcounts, taxid_Accessions, lineage):
	lineage_Readcounts = list()
	taxid_Readcounts = dict()
	for i in range(0, len(accession_Readcounts)):
		for key in taxid_Accessions:
			if accession_Readcounts[i][0] in taxid_Accessions[key]:
				if key in taxid_Readcounts:
					taxid_Readcounts[key] += accession_Readcounts[i][1]
				else:
					taxid_Readcounts[key] = accession_Readcounts[i][1]

	for key in taxid_Readcounts:
		lineage_Readcounts.append([lineage[key], taxid_Readcounts[key]])

	return lineage_Readcounts



def reverseSublists(list_of_lists):

    return [x[::-1] for x in list_of_lists]



def writeResults(outDir, lineage_Readcounts):
	results = list()
	reversed_List = list()
	resultsFile = open(outDir + '/results.tab', 'w')
	resultsFile.write('\t'.join(['Superkingdom', 'Phylum', 'Class', 'Order'\
		, 'Family', 'Genus', 'Species', 'Scientific_Name', 'TaxID', 'NumCounts']))
	resultsFile.write('\n\n')
	for sublist in range(0, len(lineage_Readcounts)):
		temp = list()
		counter = 0
		for element in lineage_Readcounts[sublist][0]:
			if element == 'NA' and counter == 0:
				temp.append(lineage_Readcounts[sublist][0][1].split('|')[1])
			elif element == 'NA' and counter > 0:
				temp.append(element)
			else:
				temp.append(element.split('|')[1] + '|' + element.split('|')[2])
			counter += 1
		if lineage_Readcounts[sublist][0][0] != 'NA':
			kept_taxid = lineage_Readcounts[sublist][0][0].split('|')[2]
		else:
			kept_taxid = lineage_Readcounts[sublist][0][1].split('|')[2]
		temp.append(kept_taxid)
		temp.append(lineage_Readcounts[sublist][1])
		results.append(temp)
	lineage_Readcounts = copy.deepcopy(results)
	for i in range(0, len(lineage_Readcounts)):
		if len(lineage_Readcounts[i]) == 9:
			lineage_Readcounts[i].insert(0, lineage_Readcounts[i][0].split('|')[0] \
				+ '|' + lineage_Readcounts[i][0].split('|')[1])
	reversed_List = reverseSublists(lineage_Readcounts)
	results_two_write = list()
	for i in range(0, len(reversed_List)):
		splitted = list()
		for j in range(1, len(reversed_List[i])):
			splitted.append(reversed_List[i][j].split('|')[0])
		splitted.insert(0, reversed_List[i][0])
		results_two_write.append(splitted)
	reversed_List = copy.deepcopy(results_two_write)
	for i in reversed_List:
		i.append(i[1])
		i.append(i[0])
		resultsFile.write('\t'.join(map(str, i[2:])))
		resultsFile.write('\n')

	return None



if __name__ == '__main__':
	parser = argparse.ArgumentParser(__file__, description = "Taxonomic Ranks")
	parser.add_argument('-quantFile', '-qf', help = 'Quantification file', default = None)
	parser.add_argument('-algorithm', '-al', help = 'Algorithm', default = 'agamemnon')
	parser.add_argument('-taxonomyFdir', '-tfd', help = 'NCBI taxonomy files directory', default = None)
	parser.add_argument('-outDir', '-od', help = 'Output directory', default = None)

	args = parser.parse_args()
	
	_accession_Readcounts = list()
	_taxid_Accessions = dict()
	_notFound = list()
	_lineage = dict()
	_lineage_Readcounts = list()
	_ranks = list()

	_accession_Readcounts = quant2list(args.quantFile, args.algorithm)
	_taxid_Accessions, _notFound = tempaccessions2taxid(_accession_Readcounts, args.taxonomyFdir + '/accessionsTaxIDs.tab')
	_lineage = buildTaxonomy(args.taxonomyFdir + '/nodes.dmp', args.taxonomyFdir + '/names.dmp', _taxid_Accessions)
	_lineage_Readcounts = counts2lineage(_accession_Readcounts, _taxid_Accessions, _lineage)
	writeResults(args.outDir, _lineage_Readcounts)
