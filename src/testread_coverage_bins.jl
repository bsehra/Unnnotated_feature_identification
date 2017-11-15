#!/usr/bin/env julia

#Created 11/1/17
#Copyright: Bhupinder Sehra

#Module ReadCoverage

using BioAlignments
#using BioSequences
#using BGZFStreams
#using GenomicFeatures

#export write_bam, header_to_array

inbam = ARGS[1]
inbamindex = ARGS[2]
binsize = parse(Int64, ARGS[3])

#Methods

function header_to_array(reader)
    #=converts @SQ information into an array containing chromosome name (::String)
    and chromosome length (::Int64)=#
    harray = []
    headr = header(reader)
    for h in headr
        if contains(string(h), "@SQ")
            h = split(string(h), "\t")
            #splits into SN:refname; LN:reflength
            refnum, chrlen = String(split(h[2], ":")[2]), parse(Int, split(h[3], ":")[2])
            push!(harray, [refnum, chrlen])
            #push!(harray, [String(split(h[2], ":")[2]), parse(Int, split(h[3], ":")[2])])
        end
    end
    return harray
end

function make_new_filename(infile, finfostr)
    #=Create new basename of file with original basename and additional string=#
    fpath = dirname(infile)
    oldfbname = splitext(basename(infile))[1]
    newfbname = oldfbname*finfostr
end


#main call

reader = open(BAM.Reader, inbam, index=inbamindex)
headarray = header_to_array(reader)

#=test
println(binsize)
println(typeof(binsize))
=#

#Iterate over input .bam file
while !eof(reader)
    for chrinfo in headarray
    #refname, reflen = headarray[index][1], headarray[index][2]; println(refname, " ", reflen)
        refname, reflen = chrinfo[1], chrinfo[2]; println(refname, " ", reflen)
        finfo = "Chr" * refname * string(binsize) * "bin_coverage"
        foutbname = make_new_filename(inbam, finfo)
        newofile = joinpath(dirname(inbam), foutbname); println(newofile, typeof(binsize))

        open(newofile, "w") do f
        #Int64(round(3/3, RoundUp))
            numsteps = Int64(round(reflen/binsize, RoundUp)); println(numsteps)
        #numsteps = Int64(round(reflen/binsize)); println(numsteps, "These are the num steps")
            readcount = 0
            for i in range(1, binsize, numsteps)
                endcount = i + binsize-1; println(endcount)
                try
                    for record in eachoverlap(reader, refname, i:endcount)
                        readcount += 1
                    end
                catch e
                    continue
                end
            #including commas turns this into a column array; leave out columns for row
                tmp = [refname "$i" "$endcount" "$readcount"]
                writedlm(f, tmp, ",") #need to format string.
                readcount = 0
            end
        end
    end
end
