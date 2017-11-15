#!/usr/bin/env python3

#Created 9/1/17
#Update 10/11/17

import sys
from sys import stderr, exit
from collections import defaultdict as dd, OrderedDict as od
#from argparse import ArgumentParser, FileType
from itertools import groupby
from operator import itemgetter
from quicksect import IntervalNode, Interval, IntervalTree

#gin = sys.argv[1]

#Methods

def dict_to_itree(fdict):
  """Convert default dictionary containing lists of feature information
  (exons, junctions, genes) to interval trees."""
   itree = IntervalTree()
   ftreedict = {}
   for key in fdict.keys():
        for vals in fdict[key]:
           #print(vals)
           itree.add(vals[1], vals[2])
        ftreedict[str(key)] = itree
        #print(ftreedict[key], "this is ftreedict[key]")
        itree = IntervalTree()
   return ftreedict


def get_anno_feature_trees(gtfin)
    genedict = dd(list)
    genes = dd(list)
    trans = od()
    exondict = dd(list)
    juncdict = dd(list)
    
    #initialise with an intervaltree
    #exondict = dd(IntervalTree())
    #juncdict = dd(IntervalTree())
    #genedict = dd(IntervalTree())

    #Iterate over .gtf file; extract 'exon' and 'gene' information lines 

    with open(gtfin, 'r') as gtf:
        """Parse valid gene/exon lines from the GTF
        file into a dict by gene_id/transcript_id"""
        for line in gtf:
            if not line.strip() or "#" in line:
                continue
            
            try:
                chrm, source, feature, left, right, score, strand, frame, values = line.split('\t')
            except ValueError:
                continue
            #Convert to 0-based coordinates - leave this in Julia as 1-based
            left, right = int(left)-1, int(right)-1
            
            if left >= right:
                continue

            if feature == 'gene':
                gvalues_dict = {}
                for attr in values.split(';')[:-1]:
                    attr, val = attr.split()
                    gvalues_dict[attr] = val.strip('"')
                if 'gene_id' not in gvalues_dict:
                    continue
                genedict[str(chrm)].append((gvalues_dict['gene_id'], left, right, strand))
            
            if feature == "exon":
                gvalues_dict = {}
                for attr in values.split(';')[:-1]:
                    attr, _, val = attr.strip().partition(' ')
                    gvalues_dict[attr] = val.strip('"')
                
                if 'gene_id' not in gvalues_dict or 'transcript_id' not in gvalues_dict:
                    continue
                transcript_id = gvalues_dict['transcript_id']

                if transcript_id not in trans:
                    trans[transcript_id] = [chrm, strand, [[left, right]]]
                    #genes[gvalues_dict['gene_id']].append(transcript_id)
                else:
                    trans[transcript_id][2].append([left, right])

    #Sort exons and merge where separating introns are <=5 bps
    #print(len(trans.items()))
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]
      
   
    #don't want unique exons per transcript...? You just want unique exons over all transcripts. 
    #Don't count shared exons with different transcript IDs.
    #so iterate over values() not items(). values are associated with keys
    tmp_exons = set()
    for chrom, strand, texons in trans.values():
        for i in range(0, len(texons)):
            tmp_exons.add((chrom, strand, texons[i][0], texons[i][1]))
    #sort by chrom and then left position
    tmp_exons = sorted(list(tmp_exons), key=itemgetter(0, 2))
    
    #Derives junctions from exons; remove duplicates by set.add()
    #Again, don't count shared junctions over several transcript IDs
    junctions = set()
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junctions.add((chrom, strand, exons[i-1][1]+1, exons[i][0]-1))
    #sort by chrom and then left position
    junctions = sorted(list(junctions), key=itemgetter(0, 2))

    """
    #python has half open intervals, exons[i-1][1] won't be counted in exon interval, and exons[1][0] won't be counted in a junction...??
    #in julia will have to do +1 and -1 as both start, end are inclusive
    junctions = set()
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junctions.add((chrom, strand, exons[i-1][1], exons[i][0]))
    #sort by chrom and then left position
    junctions = sorted(list(junctions), key=itemgetter(0, 2))    
    """



    #form dicts from lists to reduce time spent iterating over lists in generating interval trees
    for ex in tmp_exons:
        exondict[ex[0]].append(ex[1:])

    for j in junctions:
        juncdict[j[0]].append(j[1:])
    

    exonitrees = dict_to_itree(exondict)
    juncitrees = dict_to_itree(juncdict)
    geneitrees = dict_to_itree(genedict)

    #if exondict and juncdict are initialsed with an IntervalTree: NEED TO TEST!!
    """
    for ex in tmp_exons:
        exondict[ex[0]].add(ex[1], ex[2])

    for j in junctions:
        juncdict[j[0]].add(j[1], j[2])

    """




    return exonitrees, juncitrees, geneitrees
    #yield tmp_exons, junctions
    #return exondict, juncdict, genedict


#Main: Test methods





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



