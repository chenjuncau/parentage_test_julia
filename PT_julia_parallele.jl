# This is fast Julia program for PT error count. it works for large SNP panel. 
lenCount = 815  # for large 50K panel. 

dataSire = readdlm("genoSire", ' ', header=false)   
dataDam = readdlm("genoDam", ' ', header=false)   
data2 = readdlm("childGeno", ' ', header=false)    

function PT_count(x, y)
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

function PT_count_one(x, y)
    count = 0
    for mm = 1:length(y[2])
        if abs(x[2][mm] - y[2][mm]) == 2
            count += 1
        end
    end
    return count
end

now()

ped_ori = readdlm("pedigree.txt", '\t', header=false)

pedFlag = Dict()
for i = 1:size(ped_ori, 1)
    pedFlag[ped_ori[i, 1]] = i
end

gIDFlag = []
for i = 1:size(dataSire, 1)
    push!(gIDFlag, string(dataSire[i, 1]))
end

for i = 1:size(dataDam, 1)
    push!(gIDFlag, string(dataDam[i, 1]))
end

parentFlag = Dict()
for i = 1:size(dataSire, 1)
    parentFlag[dataSire[i, 1]] = i
end

now()

xxSire = []
childIndex = []
for i = 1:size(data2, 1)
    if !(string(ped_ori[pedFlag[data2[i, 1]], 2]) in gIDFlag)
        childIndex = push!(childIndex, i)
    elseif PT_count_one(data2[i, :], dataSire[parentFlag[ped_ori[pedFlag[data2[i, 1]], 2]], :]) <= lenCount
        countSire = PT_count_one(data2[i, :], dataSire[parentFlag[ped_ori[pedFlag[data2[i, 1]], 2]], :])
        ttt1 = "$(data2[i, 1]) $(ped_ori[pedFlag[data2[i, 1]], 2]) $(countSire)\n"
        xxSire = vcat(xxSire, ttt1)
        continue
    else
        childIndex = push!(childIndex, i)
    end
end

now()

yySire = []
if size(childIndex, 1) == 0
    println("There is no sire error, everything is clean\n")
else
    data22 = data2[childIndex, :]
    yySire = @parallel vcat for i = 1:size(data22, 1)
        PT_count(data22[i, :], dataSire)
    end
end

now()

parentFlag = Dict()
for i = 1:size(dataDam, 1)
    parentFlag[dataDam[i, 1]] = i
end

now()

xxDam = []
childIndex = []
for i = 1:size(data2, 1)
    if !(string(ped_ori[pedFlag[data2[i, 1]], 3]) in gIDFlag)
        childIndex = push!(childIndex, i)
    elseif PT_count_one(data2[i, :], dataDam[parentFlag[ped_ori[pedFlag[data2[i, 1]], 3]], :]) <= lenCount
        countDam = PT_count_one(data2[i, :], dataDam[parentFlag[ped_ori[pedFlag[data2[i, 1]], 3]], :])
        ttt2 = "$(data2[i, 1]) $(ped_ori[pedFlag[data2[i, 1]], 3]) $(countDam)\n"
        xxDam = vcat(xxDam, ttt2)
        continue
    else
        childIndex = push!(childIndex, i)
    end
end

now()

yyDam = []
if size(childIndex, 1) == 0
    println("There is no dam error, everything is clean\n")
else
    data22 = data2[childIndex, :]
    yyDam = @parallel vcat for i = 1:size(data22, 1)
        PT_count(data22[i, :], dataDam)
    end
end

now()

of = open("conflict_count", "w")
write(of, xxSire)
write(of, yySire)
write(of, xxDam)
write(of, yyDam)
close(of)
