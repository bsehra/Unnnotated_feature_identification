#!/usr/bin/env python3

import sys
from sys import stderr, exit
from collections import defaultdict as dd, OrderedDict as od
#from argparse import ArgumentParser, FileType
from itertools import groupby
from operator import itemgetter

#gin = sys.argv[1]

"""
def gene_separator(gline):
    return gline.split("\t")[2] == 'gene'

def get_genes(gtfin):
    glist = []
    genedict = {}
    with open(gtfin, 'r') as gtf:
        for line in gtf:
            if '#' in line or not line.strip():
                continue
            try:
                chrm, source, feature, left, right, score, strand, frame, values = line.split('\t')
            except ValueError:
                continue
            if feature == 'gene':
                left, right = int(left)-1, int(right)-1    
                gvalues_dict = {}
                for attr in values.split(';')[:-1]:
                    attr, val = attr.split()
                    gvalues_dict[attr] = val.strip('"')
                if 'gene_id' not in gvalues_dict:
                    continue
                #genedict[gvalues_dict['gene_id']] = (chrm, left, right, strand)
                genedict[str(chrm)] = (gvalues_dict['gene_id'], left, right, strand)
    return genedict

"""
def get_exons_junctions(gtfin):
    genedict = dd(list)
    genes = dd(list)
    trans = od()
    exondict = dd(list)
    juncdict = dd(list)
    
    with open(gtfin, 'r') as gtf:
        """Parse valid gene/exon lines from the GTF
        file into a dict by gene_id/transcript_id"""
        for line in gtf:
            if not line.strip() or "#" in line:
                continue
            #if '#' in line:
                #line = line.split('#')[0].strip()
            try:
                #line = line.split("\t")
                chrm, source, feature, left, right, score, strand, frame, values = line.split('\t')
            except ValueError:
                continue
            left, right = int(left)-1, int(right)-1
            
            if left >= right:
                continue
            if feature == 'gene':
                #left, right = int(left)-1, int(right)-1
                gvalues_dict = {}
                for attr in values.split(';')[:-1]:
                    attr, val = attr.split()
                    gvalues_dict[attr] = val.strip('"')
                if 'gene_id' not in gvalues_dict:
                    continue
                #genedict[gvalues_dict['gene_id']] = (chrm, left, right, strand)
                genedict[str(chrm)].append((gvalues_dict['gene_id'], left, right, strand))
            
            if feature == "exon":
                gvalues_dict = {}
                for attr in values.split(';')[:-1]:
                    #print(attr.strip().partition(' '))
                    attr, _, val = attr.strip().partition(' ')
                    gvalues_dict[attr] = val.strip('"')
                
                if 'gene_id' not in gvalues_dict or 'transcript_id' not in gvalues_dict:
                    continue
                transcript_id = gvalues_dict['transcript_id']
                if transcript_id not in trans:
                    trans[transcript_id] = [chrm, strand, [[left, right]]]
                    genes[gvalues_dict['gene_id']].append(transcript_id)
                else:
                    trans[transcript_id][2].append([left, right])

    #Sort exons and merge where separating introns are <=5 bps
    print(len(trans.items()))
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]
      
    #get unique exons
    #print(len(trans), "This is the len of trans")

    #don't want unique exons per transcript...? You just want unique exons overall, and you don't want to count shared exons with different tids as unique!
    #so iterate over values() not items()
    tmp_exons = set()
    for chrom, strand, texons in trans.values():
        for i in range(0, len(texons)):
            tmp_exons.add((chrom, strand, texons[i][0], texons[i][1]))
    tmp_exons = sorted(list(tmp_exons), key=itemgetter(0, 2))
    
    for ex in tmp_exons:
        exondict[ex[0]].append(ex[1:])

    junctions = set()
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junctions.add((chrom, strand, exons[i-1][1]+1, exons[i][0]-1))
    junctions = sorted(list(junctions), key=itemgetter(0, 2))

    for j in junctions:
        juncdict[j[0]].append(j[1:])
    
    #yield tmp_exons, junctions
    return exondict, juncdict, genedict


"""
from itertools import islice

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

n_items = take(30, sortex)
print(n_items)

n_items = take(30, sortjunc)
print(n_items)

print(take(30, transcripts.items()))

"""

"""
if __name__ == '__main__':
    from sys import stderr, exit
    import sys
    from collections import defaultdict as dd, OrderedDict as od
    from operator import itemgetter
    get_exons_junctions(sys.argv[1])
    if not sys.argv[1]:
        print("no argument")
        exit(1)
"""
