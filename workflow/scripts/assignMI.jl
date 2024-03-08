#!/usr/bin/env julia

try
    using XAM ;
    using BGZFStreams ;
catch
    using Pkg; Pkg.add(["XAM", "BGZFStreams"],  io=devnull) ;
    using XAM ;
    using BGZFStreams ;
end


function getBX(itr::BAM.Record, fallback::String)::String
    try
        itr["BX"]::String
    catch
        fallback::String
    end
end

function getMI(itr::BAM.Record, fallback::Int64)::Int64
    try
        Int64(itr["MI"])::Int64
    catch
        fallback::Int64
    end
end

function write_validbx(bam, record::BAM.Record, molID::Int64)
    """
    bam: the output bam
    alnrecord: the pysam alignment record
    molID: the "mol_id" entry of a barcode dictionary
    Formats an alignment record to include the MI tag
    and the BX at the end and writes it to the output
    bam file. Replaces existing MI tag, if exists.
    """
    # get all the tags except MI b/c it's being replaced (if exists)
    # also remove DI because it's not necessary
    tags = BAM.auxdata(record)
    tags["MI"] = molID
    # find which tag index is the BX tag
    # write record to output file
    bam.write(alnrecord)
end
function write_invalidbx(bam, alnrecord)
    """
    bam: the output bam
    alnrecord: the pysam alignment record
    Formats an alignment record to include the BX 
    at the end and writes it to the output
    bam file. Removes existing MI tag, if exists.
    """
    # get all the tags except MI b/c it's being replaced (if exists)
    # this won't write a new MI, but keeping an existing one
    # may create incorrect molecule associations by chance
    # also remove DI because it's not necessary
    tags = [j for i,j in enumerate(alnrecord.get_tags()) if j[0] not in ['MI', 'DI']]
    # find which tag index is the BX tag
    BX_idx = [i for i,j in enumerate(tags) if j[0] == 'BX'][0]
    # get the list of indices for the tag list
    idx = [i for i in range(len(tags))]
    # if BX isn't already last, make sure BX is at the end
    if tags[-1][0] != 'BX'
        # swap it with whatever is last
        tags[BX_idx], tags[-1] = tags[-1], tags[BX_idx]
        # update the record's tags
        alnrecord.set_tags(tags)
    end
    # write record to output file
    bam.write(alnrecord)
end
function write_missingbx(bam, alnrecord)
    """
    bam: the output bam
    alnrecord: the pysam alignment record
    Formats an alignment record to include invalid BX 
    at the end and writes it to the output
    bam file. Removes existing MI tag, if exists.
    """
    # get all the tags except MI b/c it's being replaced (if exists)
    # this won't write a new MI, but keeping an existing one
    # may create incorrect molecule associations by chance
    # also remove DI because it's not necessary
    # removes BX... just in case. It's not supposed to be there to begin with
    tags = [j for i,j in enumerate(alnrecord.get_tags()) if j[0] not in ['MI', 'DI', 'BX']]
    tags.append(("BX", "A00C00B00D00"))
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)
end
    

bam_input = snakemake.input[0]
bamw = BAM.Writer(BGZFStream(open("my-data.bam", "w"), "w"))
# initialize the dict
d = dict()
# chromlast keeps track of the last chromosome so we can
# clear the dict when it's a new contig/chromosome
chromlast = False
# MI is the name of the current molecule, starting a 1 (0+1)
MI = 0
write(bamw, rec) # rec = BAM.Record
if os.path.exists(bam_input) and bam_input.lower().endswith(".sam"):
    alnfile = pysam.AlignmentFile(bam_input)
elif os.path.exists(bam_input) and bam_input.lower().endswith(".bam"):
    if os.path.exists(bam_input + ".bai"):
        alnfile = pysam.AlignmentFile(bam_input)
    else:
        print(f"Error: {bam_input} requires a matching {bam_input}.bai index file, but one wasn\'t found.", file = sys.stderr)
        exit(1)
else:
    print(f"Error: {bam_input} not found", file = sys.stderr)
    exit(1)

# iniitalize output file
#alnfile = pysam.AlignmentFile("/home/pdimens/Documents/harpy/test/bam/sample1.bam")
outfile = pysam.AlignmentFile(snakemake.output[0], "wb", template = alnfile)
#outfile = pysam.AlignmentFile("/home/pdimens/Documents/harpy/test/bam/test.bam", "w", template = alnfile)

for record in alnfile.fetch():
    chrm = record.reference_name
    bp   = record.query_alignment_length
    # check if the current chromosome is different from the previous one
    # if so, empty the dict (a consideration for RAM usage)
    if chromlast != False and chrm != chromlast:
        d = dict()
    if record.is_unmapped:
        # skip, don't output
        chromlast = chrm
        continue

    try:
        bx = record.get_tag("BX")
        # do a regex search to find X00 pattern in the BX
        if re.search("[ABCD]0{2,4}", bx):
            # if found, invalid
            write_invalidbx(outfile, record)
            chromlast = chrm
            continue
    except:
        # There is no bx tag
        write_missingbx(outfile, record)
        chromlast = chrm
        continue
    
    aln = record.get_blocks()
    if not aln:
        # unaligned, skip and don't output
        chromlast = chrm
        continue

    # logic to accommodate split records 
    # start position of first alignment
    pos_start = aln[0][0]
    # end position of last alignment
    pos_end   = aln[-1][1]

    # create bx entry if it's not present
    if bx not in d.keys():
        # increment MI b/c it's a new molecule
        MI += 1 
        d[bx] = {
            "lastpos" : pos_end,
            "current_suffix": 0,
            "mol_id": MI
        }
        # write and move on
        write_validbx(outfile, record, d[bx]["mol_id"])
        chromlast = chrm
        continue

    # store the original barcode as `orig` b/c we might need to suffix it
    orig = bx
    # if there is a suffix, append it to the barcode name
    if d[orig]["current_suffix"] > 0:
        bx = orig + "." + str(d[orig]["current_suffix"])

    # distance from last alignment = current aln start - previous aln end
    dist = pos_start - d[bx]["lastpos"]
    # if the distance between alignments is > cutoff, it's a different molecule
    # so we'll +1 the suffix of the original barcode and relabel this one as 
    # BX + suffix. Since it's a new entry, we initialize it and move on
    if dist > snakemake.params[0]:
        # increment MI b/c it's a new molecule
        MI += 1 
        # increment original barcode's suffix
        d[orig]["current_suffix"] += 1
        bx = orig + "." + str(d[orig]["current_suffix"])
        # add new entry for new suffixed barcode with unique MI
        d[bx] = {
            "lastpos" : pos_end,
            "current_suffix": 0,
            "mol_id": MI
        }
        # write and move on
        write_validbx(outfile, record, d[bx]["mol_id"])
        chromlast = chrm
        continue 

    if record.is_reverse or (record.is_forward and not record.is_paired):
        # set the last position to be the end of current alignment
        d[bx]["lastpos"] = pos_end

    # if it hasn't moved on by now, it's a record for an
    # existing barcode/molecule. Write the record.
    write_validbx(outfile, record, d[bx]["mol_id"])

    # update the chromosome tracker
    chromlast = chrm

alnfile.close()
outfile.close()

# index the output file
pysam.index(snakemake.output[0])

# exit gracefully
exit(0)