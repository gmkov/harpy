#! /usr/bin/env julia

try
    using XAM
catch
    using Pkg; Pkg.add("XAM")
    using XAM
end

function parseBX!(rec, inv_dict, rgx)
    bx = rec["BX"]::String
    barcode = match(rgx, bx).match
    let invalid::Int64 = 0
        for i in ["A","B","C","D"]
            check = Int64(occursin(i * "00", barcode))
            inv_dict[i] += check
            invalid += check
        end
        return iszero(invalid)
    end
end

# setup
haplotag_regex = r"A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}"
reader = open(BAM.Reader, ARGS[1])
record = BAM.Record()

let n_reads = 0, n_bx = 0, n_valid = 0, invalid_dict = Dict{String, Int64}("A" => 0, "B" => 0, "C" => 0, "D" => 0)
    while !eof(reader)
        empty!(record)
        read!(reader, record)
        n_reads += 1
        try
            valid = parseBX!(record, invalid_dict, haplotag_regex)
            n_bx += 1
            n_valid += valid
        catch
            continue
        end
    end
    open("test.txt", "w") do fout
        write(fout, "totalReads\t" *"$n_reads" * "\n")
        write(fout, "bxTagCount\t" * "$n_bx" * "\n")
        write(fout, "bxValid\t" * "$n_valid" * "\n")
        write(fout, "bxInvalid\t" * string(n_bx - n_valid) * "\n")
        write(fout, "A00\t" * string(invalid_dict["A"])* "\n")
        write(fout, "C00\t" * string(invalid_dict["C"])* "\n")
        write(fout, "B00\t" * string(invalid_dict["B"])* "\n")
        write(fout, "D00\t" * string(invalid_dict["D"])* "\n")
    end
end