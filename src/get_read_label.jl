#!/usr/bin/env python

#Created 9/29/17
#Update 11/4/17
#Bhupinder Sehra

#=Information on the alignment of each read from a given .bam/.sam (representing 1 biorep) is provided in an input file 'readinfo'
created using a Julia script (v0.6). The script iterates over each read and provides a label identifying whether the alignment of each
read is 'consistent' with given annotation (True) or 'inconsistent' (False).

Inconsistent reads are potentially indicative of the presence of unannotated features, including novel splice junctions/exon transcriptional 
boundaries or other genomic phenomena.=#

using BioAlignments
using BGZFStreams
using GenomicFeatures

include("./get_annotated_features.jl")
using get_annotated_features
#using ResumableFunctions - incorporate to use python/C# style macros


gfin = ARGS[1]
inbam = ARGS[2]
inbamindex = ARGS[3]
#to write out
classfile = ARGS[4]


etreedict, jtreedict, gtreedict = get_anno_feature_trees(gfin)

#julia version? BAM.cigar_rle(record): print out see what this looks like
function get_gap_coord(cigtup, refst)
  #=Find genomic coordinates of gap ('N') in a read. 
  cigtup is BAM.cigar_rle(record): [(ops), (opslen)].
  Tested - works!=#
  ops = cigtup[1]
  oplen = cigtup[2]
  rst, rend = 0, 0
  incops = "MDX="
  gapcoords = Array{Int64, 1}()

  for (index, value) in enumerate(ops)
    #boolean test; value = operation
    value = string(value) #converting BioAlignments.Operation -> string() works.
    if contains(incops, value)
        rst += Int64(oplen[index])
    elseif value == "N"
        rst = rst + refst - 1 
        rend  = rst + Int64(oplen[index])
        push!(gapcoords, (rst, rend))    
        rst = rend - refst + 1
    end
  end
return gapcoords
end

#merge two methods check_njintervals and check_genes 
function check_njinterval(rst, rend, alen, ftree, **checkgenes):
    #=Check if an ungapped read is consistent with annotated features; 
    interval checking of read coordinates with exon/gene interval tree.
    Consistent -> T; inconsistent -> F
    keyword (boolean) checkgenes=#

    #returns an array of intervals
    olap_label = nothing
    inter = collect(intersect(ftree, (rst, rend))) 

    #if read overlaps with only one exon across full length of read -> true
    if length(inter) == 1
        ist, iend = inter[1].first, inter[1].last
        olap = minimum([iend-ist, iend-rst, rend-rst, rend-ist])
       #Evaluate by if check_njinterval() - True if it is consistent. If no check_njinternal - is of interest
        olap_label = alen-2 <= olap <= alen+2

    #Short reads: max overlap should be length of read; long reads, max overlap is length of an exon interval
    #If max overlap != alignment length -> true
    elseif length(inter) > 1
        #list comprehension over all intervals; find max overlap
        #max_olap = maximum([minimum(inter[i].last-inter[i].first, inter[i].last-rst, rend-rst, rend-inter[i].first) for i in 1:length(inter)]) - from python
        max_olap = maximum(minimum([inter[i].last-inter[i].first, inter[i].last-rst, rend-rst, rend-inter[i].first]) for i in 1:length(inter)) #julia syntax
        olap_label = alen-2 <= max_olap <= alen+2

    #if read not found to overlap with any features
    elseif length(inter) == 0
        if checkgenes == false
          olap_label = "None" #checking exons
        else:
          olap_label = true #result of checking genes
        end
    end
return olap_label
end    

#g = ( (a, b) for a in v1, b in v2 if a+b > 7 )
#maximum(minimum([inter[i].last-inter[i].first, inter[i].last-rst, rend-rst, rend-inter[i].first]) for i in 1:length(inter))
#minimum(minimum([1+i 2+i 3+i 4+i]) for i in 1:10) - this seems to work



function check_jinterval(gapst, gapend, ftree):
    #=interval checking given an input interval tree: st and end are the coordinates of non-gapped read=#
    #obool = nothing
    inter = collect(intersect(ftree, (gapst, gapend)))
    
    if length(inter) == 1:
        ist, iend = inter[1].first, inter[1].last      
        obool = (ist-2 <= gapst <= ist+2) && (iend-2 <= gapend <= iend+2)

    #check each interval for consistency with annotated junction coordinates
    #if any intervals are consistent with gap coordinates -> true
    elseif length(inter) > 1:  
        obool = any((inter[i].first-2 <= gapst <= inter[i].first+2 && inter[i].last-2 <= gapend <= inter[i].last+2) for i in 1:length(inter)))

    #gap alignment does not overlap with any annotated junctions
    elseif length(inter) == 0
        obool = false #no junction meets annotation
    end
return obool
end


function check_jinterval(gapst, gapend, ftree):
    #=interval checking given an input interval tree: st and end are the coordinates of non-gapped read=#
    #try/catch error: if gap in spliced read maps to a region with no annotated junctions -> false (not consistent with annotation)
    #check each interval for consistency with annotated junction coordinates
    #if any intervals are consistent with gap coordinates -> true
    #obool = nothing
    inter = collect(intersect(ftree, (gapst, gapend)))
    if length(inter) == 1:
        ist, iend = inter[1].first, inter[1].last      
        obool = (ist-2 <= gapst <= ist+2) && (iend-2 <= gapend <= iend+2)
    
    elseif length(inter) > 1:  
        obool = any(inter[i].first-2 <= gapst <= inter[i].first+2 && inter[i].last-2 <= gapend <= inter[i].last+2 for i in 1:length(inter)))

    #gap alignment does not overlap with any annotated junctions
    elseif length(inter) == 0
        obool = false #no junction meets annotation
    end
return obool
end


#try using ResumableFunctions: @yield macro to use generator instead of return 

function get_class(cigar, cigtuples, rdst, rdend, alignlen, refnum, ftree, gtreedict):
  #if read is ungapped, check exon interval tree; if no overlap found ("None") check gene intervals
  #else yield checkint (consistent -> true; inconsistent -> false)
  #rlabel = nothing
  if contains(cigar, "N")
      #@yield rlabel
      #checkint = check_njinterval(rdst, rdend, alignlen, exontree)
      #if checkint == "None": #result of no exons overlapped - check genes and yield result
      if check_njinterval(rdst, rdend, alignlen, exontree) == "None"
          gtree = gtreedict[refnum] #only derive gtree when needed 
          rlabel = check_njinterval(rdst, rdend, alignlen, gtree, checkgenes=True)
      else:
          rlabel = checkint
      end
  #process gapped reads; method accounts for multiple gaps; if any gaps are inconsistent with annotation -> false
  elseif !contains(cigar, "N")
      gapc = get_gap_coord(cigtuples, rdst)
      rlabel = all(check_jinterval(gap[0], gap[1], ftree) for gap in gapc) # -> true only if all are true
      #any((5 <= i <= 8 && 2 <= i <= 5) == true for i in 1:10) - this works!
  end
#rlabel
return rlabel
end

#Create dictionaries containing annotated exon, junction and gene information from input GTF file
#change the original module to go straight to producing interval trees

#Main program

#Iterate over .bam/.sam file(s); generate training data set
#first file (Single end reads) - input array of bam files; labelfile in append mode

#bamarray = [inbam1, inbam2, inbam3]
#for bamf in bamarray:


#This is memory efficient but slow! Increase number of threads
reader = open(BAM.Reader, inbam, index=inbamindex)
#=Take in mapped reads from a bam file and write out select information to file; adding cigar_rle tuple=#
#open(labelfile, "a") do lf #or append mode???
open(labelfile, "w") do lf #or append mode???
  record = BAM.Record()
  while !eof(reader)
    read!(reader, record)
    if BAM.ismapped(record)
          #tmp = [BAM.tempname(record) atype biorep BAM.mappingquality(record) BAM.flag(record) BAM.cigar(record) BAM.position(record) BAM.rightposition(record) BAM.alignlength(record) BAM.refname(record)]
        tmparray = [BAM.tempname(record), biorep, BAM.mappingquality(record), BAM.flag(record), BAM.cigar(record), BAM.cigar_rle(record), BAM.position(record), BAM.rightposition(record), BAM.alignlength(record), BAM.refname(record)]              
        if !contains(cigar, "N")
            featuretree = extreed[string(BAM.refname)]                      
            #genetree = gtreed[string(BAM.refname)]
        else
            #some chromosomes have only 1-exon transcripts and no junctions 
            try
                featuretree  = jtreed[string(BAM.refname)] 
            catch KeyError
                #initialise with an empty tree
                featuretree =  IntervalTree{Int, Interval{Int}}()
            end
        end
        classbool = get_class(BAM.cigar(record), BAM.cigar_rle(record), BAM.position(record), BAM.rightposition(record), BAM.alignlength(record), BAM.refname(record), featuretree, gtreedict)
        append!(tmp, classbool)
        if lfiletype == "bed"
            description = string(BAM.tempname(record)) * "_" * string(cigar) * "_" * string(Int8(classbool)) #how to concatenate string??
            bedarray = [BAM.refname(record) BAM.position(record), BAM.rightposition(record), description]
            writedlm(lf, bedarray, "\t")
        else
            writedlm(lf, tmparray, "\t")
    end
  end  
end

#=
open(readinfo, "r") do rinf
  open(classfile, "w") do outf:
    if bed == true
        write(outf, track + " " + description + "\n")
    end

    for line in rinf:
      if not line or not line.strip():
          continue
      end
    
      line = strip(line, " ")
        #extract information to get label
      rid, atype, biorep, mapq, flag, cigar, stcoord, endcoord, alen, ref = line.split("\t")
      stcoord, endcoord, alen, refnum = Int64(stcoord), Int64(endcoord), Int64(alen), string(ref)

    #Extract interval tree for relevant chromosome. Exception handling if key (chromosome/refnum is missing missing)
      try:
          extree, genetree, junctree = etreedict[ref], jtreedict[ref], gtreedict[ref]
      catch e:
          pass
      end

    #processing each read: get label (T/F)
      rbool = Bool(get_class(cigar, stcoord, endcoord, alen, extree, junctree, genetree)) 
   
    #write out results to file (same layout as input file with label as extra column or .bed file)
      if not bed:
          write(outf, line + "\t" + string(Int8(rbool)) + "\n")
      else:
        #write BED file
        if rbool == true
            write(outf, refnum + "\t" + str(start) + "\t" + str(end) + "\t" + str(rid) + "_" + str(int(rbool)) + cigar + "\n")         
        end
        #counter; keep track of lines
        i += 1
        println(i)
      end
  end
end
print("done!")
=#































