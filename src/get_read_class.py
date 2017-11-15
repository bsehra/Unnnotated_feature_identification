#!/usr/bin/env python3

#Created 9/29/17
#Update 11/4/17
#Bhupinder Sehra

"""Information on the alignment of each read from a given .bam/.sam (representing 1 biorep) is provided in an input file 'readinfo'
created using a Julia script (v0.6). The script iterates over each read and provides a label identifying whether the alignment of each
read is 'consistent' with given annotation (True) or 'inconsistent' (False).

Inconsistent reads are potentially indicative of the presence of unannotated features, including novel splice junctions/exon transcriptional 
boundaries or other genomic phenomena.""" 


import sys
from Getgtfinfotest import get_exons_junctions
from collections import defaultdict as dd
from quicksect import IntervalNode, Interval, IntervalTree
import re

gfin = sys.argv[1] #use Araport11
readinfo = sys.argv[2]
classfile = sys.argv[3]

def dict_to_itree(fdict):
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


#julia version? BAM.cigar_rle(record)
def cigarstr_tuples(cigsc):
    """Output a tuple rep of cigar score."""
    return re.findall('(\d+)([MIDNSHPX])?', cigsc)


#julia version
#need to allow for 'wiggle room of ~1-2bp or so?'
def get_gap_coord(cig_tup, refst):
    """Find coordinates of the start and end of a gap within a splice read given the left genomic coordinate."""
    st, end = 0, 0
    incops = "MDX="
    stopops = "ISHP"
    gapcoords = []
    for tup in cig_tup:
        #if M,D,X or += advance count
        if tup[1] in incops:
            st += int(tup[0])
        #ignore deletions in the reference, soft/hard clipping
        if tup[1] in stopops:
            pass
        if tup[1] == "N":
            st = refst + st - 1
            end = st + int(tup[0]) - 1
            gapcoords.append((st, end))
            st = end - refst + 1
    return gapcoords


"""True = consistent; False = inconsisten"""


#merge two methods check_njintervals and check_genes 
def check_njinterval(st, end, alen, ftree, **checkgenes):
    """Check if an ungapped read is consistent with annotated features; 
    interval checking of read coordinates with exon/gene interval tree.
    Consistent -> T; inconsistent -> F
    keyword (boolean) checkgenes"""

    #inter = IntervalTree()
    inter = ftree.search(st, end)
    #if read overlaps with only one exon across full length of read -> True
    if len(inter) == 1:
       ist, iend = inter[0].start, inter[0].end
       olap = min(iend-ist, iend-st, end-st, end-ist)
       #Evaluate by if check_njinterval() - True if it is consistent. If no check_njinternal - is of interest
       return alen-2 <= olap <= alen+2

    #Short reads: max overlap should be length of read; long reads, max overlap is length of an exon interval
    #If max overlap != alignment length -> True
    if len(inter) > 1:
        max_olap = max([min(inter[i].end-inter[i].start, inter[i].end-st, end-st, end-inter[i].start) for i in range(len(inter))])
        return alen-2 <= max_olap <= alen+2

    #if read not found to overlap with any features
    if not inter:
      if not checkgenes:
        return "None" #checking exons
      else:
        return True #result of checking genes


#rewrite this as if ist-2 <= st <= ist+2 and iend-2 <= end <= iend+2 -> False
       #or return ist-2 => st >= ist+2 and iend-2 => end >= iend+2 ->  True : do this in julia script
def check_jinterval(st, end, ftree):
    """interval checking given an input interval tree: st and end are the coordinates of non-gapped read"""
    #inter = IntervalTree()
    inter = ftree.search(st, end)
    if len(inter) == 1:
      ist, iend = inter[0].start, inter[0].end       
      return (ist-2 <= st <= ist+2) and (iend-2 <= end <= iend+2)
    
    #check each interval for consistency with annotated junction coordinates
    #if any intervals are consistent with gap coordinates, then return True
    if len(inter) > 1:  
      return any((inter[i].start-2 <= st <= inter[i].start+2) and (inter[i].end-2 <= end <= inter[i].end+2) for i in range(len(inter)))

    #gap alignment does not overlap with any annotated junctions
    if not inter:
       return True


#This may be redundant as a method; combine with check_njinterval?
def check_genes(st, end, alen, gtree):
    """Given an interval tree of gene coordinates and genomic coordinates of an aligned read, determine whether
     the read falls completely within any gene interval (F) or does not (partial or no overlap; T)."""
    genes = gtree.search(st, end)
    if len(genes) == 1:
      ist, iend = inter[0].start, inter[0].end   
      olap = min(iend-ist, iend-st, end-st, end-ist)
      return alen-2 <= olap <= alen+2
      
    #Falls within several gene boundaries; if read is aligned fully within any gene boundary (intervals) -> False 
    if len(genes) > 1:
       max_olap = max([min(genes[i].end-genes[i].start, genes[i].end-st, end-st, end-genes[i].start) for i in range(len(genes))])
      return alen-2 <= max_olap <= alen+2
       
    #does not fall within any annotated gene boundary
    if not genes:
      return True


def get_class(cigar, rdst, rdend, alignlen, etree, jtree, gtree):
  #if not in("N", cigar) - in julia
  #from main method
  #if read is ungapped, check exon interval tree; if no overlap found ("None") check gene intervals
  #else yield checkint (consistent -> true; inconsistent -> False)
  if "N" not in cigar:
      checkint = check_njinterval(rdst, rdend, alignlen, etree)
      if checkint == "None": #result of no exons overlapped - check genes and yield result
          yield check_njinterval(rdst, rdend, alignlen, gtree, checkgenes=True)
      else:
          yield checkint  

  #from main method below
  #process gapped reads; method accounts for multiple gaps
  if "N" in cigar:
      rbool = None
      cigtup = cigarstr_tuples(cigar)
      gapc = get_gap_coord(cigtup, rdst)
      yield any(check_jinterval(gap[0], gap[1], jtree) for gap in gapc)

  
#Create dictionaries containing annotated exon, junction and gene information from input GTF file
#change the original module to go straight to producing interval trees
exondict, juncdict, genedict = get_exons_junctions(gfin)

#convert features (exons, junctions, genes) specific to each chromosome to interval trees; store in dictionary
exonitrees = dict_to_itree(exondict)
juncitrees = dict_to_itree(juncdict)
geneitrees = dict_to_itree(genedict)

#clear unwanted dictionaries
exondict.clear()
juncdict.clear()
genedict.clear()


#Main program

#variables
i = 0
rbool = None
bed = True
track = "track name=" + '"Mouse_STAR_SE_True"'
description = "description=" + '"Truereads"'

#Iterate over file containing key information extracted for each read from .bam/.sam file
#File was produced in Julia
with open(readinfo, "r") as rinf, open(classfile, "w") as outf:
  if bed:
    outf.write(track + " " + description + "\n")

  for line in rinf:
    if not line or not line.strip():
      continue
    line = line.strip()
        #extract information to get label
    rid, atype, biorep, mapq, flag, cigar, start, end, alen, refnum = line.split("\t")
    start, end, alen, refnum = int(start), int(end), int(alen), str(refnum)

    #Extract interval tree for relevant chromosome. Exception handling if key (chromosome/refnum is missing missing)
    try:
        extree, genetree, junctree = exonitrees[refnum], geneitrees[refnum], juncitrees[refnum]
    except KeyError:
      pass

    #processing each read: get label (T/F)
    rbool = bool(get_class(cigar, start, end, alen, extree, junctree, genetree)) 
   
    #write out results to file (same layout as input file with label as extra column or .bed file)
    if not bed:
        outf.write(line + "\t" + str(int(rbool)) + "\n")
    else:
        #write BED file
        if rbool:
            outf.write(str(refnum) + "\t" + str(start) + "\t" + str(end) + "\t" + str(rid) + "_" + str(int(rbool)) + cigar + "\n")         
        #counter; keep track of lines
        i += 1
        print(i)
print("done!")





























#order per line in readinfo file:
#[BAM.tempname(record) atype biorep BAM.mappingquality(record) BAM.flag(record) BAM.cigar(record) BAM.position(record) BAM.rightposition(record) BAM.alignlength(record) BAM.refname(record)]
"""                   
with open(readinfo, "r") as rinf, open(classfile, "w") as outf:
    for line in rinf:
        if not line or not line.strip():
            continue
        line = line.strip()
        #extract information to get label
        rid, atype, biorep, mapq, flag, cigar, start, end, alen, refnum = line.split("\t")
        start, end, alen, refnum = int(start), int(end), int(alen), str(refnum)
        try:
           extree, genetree, junctree = exonitrees[refnum], geneitrees[refnum], juncitrees[refnum]
        except KeyError:
           pass
        
        rbool = bool(get_class(cigar, start, end, alen, extree, junctree, genetree))
        #print(rbool, "this is bool")
        outf.write(line + "\t" + str(int(rbool)) + "\n")
        #print(line, " ", rbool)
print("done!")
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





















