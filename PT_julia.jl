# PT program using Julia. 
using Distributed

# Read genotype data
data1 = readdlm("parentGeno", ' ', header=false)
data2 = readdlm("childGeno", ' ', header=false)

# Define function for counting differences
@everywhere function PT_count(x, y)
    z = []
    for j = 1:size(y, 1)
        count = 0
        for mm = 1:length(y[1, 2])
            if abs(x[2][mm] - y[j, 2][mm]) == 2
                count += 1
            end
        end
        t = "$(x[1]) $(y[j, 1]) $(count)\n"
        z = vcat(z, t)
    end
    return z
end

# Parallel computation of PT_count
yy = @distributed vcat for i = 1:size(data2, 1)
    PT_count(data2[i, :], data1)
end

# Write results to file
open("count_of.txt", "w") do of
    write(of, yy)
end

# Parse parent gender data
gID = Dict()
open("parent_gender") do f
    for ln in eachline(f)
        blupID, sex = split(ln, r" ")
        gID[blupID] = chomp(sex)
    end
end

# Parse parent genotype data
gIDFlag = []
open("parentGeno") do f
    for ln in eachline(f)
        ID1, geno = split(ln, r" ")
        push!(gIDFlag, string(ID1))
    end
end

# Read and process pedigree data
ped_ori = readdlm("ped_database.txt", ' ', header=false)
genoPed = ped_ori[:, 1:6]
genoPed[:, 2:6] .= 0

genoPedIndex = Dict()
for i = 1:size(genoPed, 1)
    genoPedIndex[genoPed[i, 1]] = i
end

open("count_of.txt") do f
    for ln in eachline(f)
        ID1, ID2, errorC = split(ln, " ")
        if parse(Int64, chomp(errorC)) > 20  # 2K SNP chip
            continue
        end
        if typeof(ID2) != String
            ID2 = string(ID2)
        end
        if typeof(ID1) != String
            ID1 = string(ID1)
        end
        if gID[ID2] == "M"
            genoPed[genoPedIndex[ID1], 2] = ID2
            genoPed[genoPedIndex[ID1], 4] += 1
            genoPed[genoPedIndex[ID1], 6] += 1
        end
        if gID[ID2] == "F"
            genoPed[genoPedIndex[ID1], 3] = ID2
            genoPed[genoPedIndex[ID1], 5] += 1
            genoPed[genoPedIndex[ID1], 6] += 1
        end
    end
end

# Write geno pedigree file
open("geno_ped.txt", "w") do of
    for i = 1:size(genoPed, 1)
        write(of, "$(genoPed[i, 1]) $(genoPed[i, 2]) $(genoPed[i, 3]) $(genoPed[i, 4]) $(genoPed[i, 5]) $(genoPed[i, 6])\n")
    end
end
