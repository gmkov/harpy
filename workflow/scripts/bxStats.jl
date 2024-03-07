#! /usr/bin/env julia

try
    using XAM
    using GZip
catch
    using Pkg; Pkg.add(["XAM", "GZip"])
    using XAM
    using GZip
end

## worker functions ##
function newInvalidDict()::Dict{String, Int64}
    Dict{String, Int64}("start" =>  0, "end"=> 0, "bp" => 0, "n" => 0, "lastpos" => 0,"mindist" => 0)
end

function newValidDict(startpos::Int64, endpos::Int64, bp::Int64)::Dict{String, Int64}
    Dict{String, Int64}("start" =>  startpos, "end" => endpos, "bp" => bp, "n" => 1, "lastpos" => endpos, "mindist" => -1)
end

function getBX(itr::BAM.Record, fallback::String)::String
    try
        itr["BX"]::String
    catch
        fallback::String
    end
end

function getMI(itr::BAM.Record, fallback::String)::Int64
    try
        itr["MI"]::Int64
    catch
        fallback::Int64
    end
end

isforward(record::BAM.Record)::Bool = BAM.flags(record) & 0x10 == 0
ispaired(record::BAM.Record)::Bool = BAM.flags(record) & 0x80 == 0
isduplicate(record::BAM.Record)::Bool = BAM.flags(record) & 0x400 == 0

function updatedict!(d::Dict{String, Int64}, record::BAM.Record)::Bool
    """
    Safely pull the BX and MI tags from the alignment record and either
    update the dictionary entry for the MI tag or create a new one.
    """
    bp::Int64  = BAM.alignlength(record)
    bx::String = getBX(record, "noBX")
    mi::Int64  = getMI(record, -1)
    if bx != "noBX"
        bad_beadtag = occursin(r"[ABCD]0{2,4}", bx)
        bx = bad_beadtag ? "invalidBX" : bx
    end
    valid = !occursin(r"[ABCD]0{2,4}", bx) && bx != "noBX"
    if !valid
        get!(d, mi, newInvalidDict())
        d[mi]["bp"] += bp
        d[mi]["n"] += 1
        return false
    end
    pos_start = BAM.position(record)
    pos_end = BAM.rightposition(record) 
    if mi ∉ keys(d)
        d[mi] = newValidDict(pos_start, pos_end, bp)
        return false         
    end

    # only calculate the minimum distance between alignments
    # if it's a forward read or an unpaired reverse read
    dist = pos_start - d[mi]["lastpos"]
    if isforward(record) || (!ispaired(record) && !isforward(record))
        _dist = d[mi]["mindist"]
        d[mi]["mindist"] = (dist < _dist || _dist < 0) ? dist : _dist

    # update the basic alignment info of the barcode
    d[mi]["bp"] += bp
    d[mi]["n"]  += 1

    # only if low < currentlow or high > currenthigh
    if pos_start < d[mi]["start"]
        d[mi]["start"] = pos_start
    end
    if pos_end > d[mi]["end"]
        d[mi]["end"] = pos_end
    end
    if !isforward(record) || (isforward(record) && !ispaired(record))
        # set the last position to be the end of current alignment
        d[mi]["lastpos"] = pos_end
    end
    return false
end

# define write function
# it will only be called when the current alignment's chromosome doesn't
# match the chromosome from the previous alignment
function writestats(x::T,chr::String) where T<:Dict
    for (mi,stats) in x
        stats["inferred"] = stats["end"] - stats["start"] 
        stats["mindist"] = max(0, stats["mindist"])
        outtext = "$chr\t$mi\t" * join([stats[i] for i in ["n", "start","end", "inferred", "bp", "mindist"]], "\t")
        write(outfile, outtext * "\n")
    end
end

write(outfile, "contig\tmolecule\treads\tstart\tend\tlength_inferred\taligned_bp\tmindist\n")


####
# pick an unreasonable string you shouldnt name contigs after (Unicode \0)
# mi = -1 if missing
let chromlast = "\u0", d = Dict{Int64,Dict{String,Int64}}()
    reader = open(BAM.Reader, ARGS[1])
    outfile = GZip.open(ARGS[2], "w")
    record = BAM.Record()
    while !eof(reader)
        empty!(record)
        read!(reader, record)
        chrm::String = BAM.refname(record)
        bp::Int64    = BAM.alignlength(record)
        # check if the current chromosome is different from the previous one
        # if so, print the dict to file and empty it (a consideration for RAM usage)
        if chrm != chromlast && chromlast != "\u0"
            writestats(d, chromlast)
            empty!(d)
        end
        chromlast = chrm
        if isduplicate(record) || !ismapped(record)
            continue
        end
        # create up update an MI dict entry
        updatedict!(d, record)
    end
end