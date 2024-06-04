# Jun Chen
# 11-21-2017
# this one is good for HV41
# PT program to find the parents. 
# check all the chicken.
# very nice.
# /biotech/genomic_selection/HV-58/run/mg222/PT_2017_11_18_01_01    # test here. 
# chenju@cvirdlux04> grep ^222 geno4PT > k1
# chenju@cvirdlux04> grep -v ^222 geno4PT > k2
 
# Z:\project_gs_pipe\source\convRawGenoAB-012.jl # using this one. 
awk '{print $1,0,0,0,0,0}' childGeno > ped_database.txt

# grep -v ^243 geno4PT > parentGeno
# cd /biotech/genomic_selection/DC-88/run/mg58/PT_2008030807
# grep -v ^DC88_58 geno4PT > parentGeno
# grep ^DC88_58 geno4PT > childGeno


julia -p 30     #  30 cpus only 6 mins, 15 cpus only 15 mins.   # very good. 
# need parentGeno,childGeno,parent_gender,ped_database.txt (original pedigree)
# ---------------method 1.  parallele------------------------------------ 
# cd("/biotech/genomic_selection/HV-58/run/mg222/PT_2017_11_18_01_01")

data1=readdlm("parentGeno",' ',header=false)   # k2 
data2=readdlm("childGeno",' ',header=false)    # This is for the curent MG.  k1

# data1=readdlm("no228",' ',header=false)   # k2 
# data2=readdlm("232",' ',header=false)    # This is for the curent MG.  k1

# addprocs(10)
@everywhere function PT_count(x,y)     # it is work in V5, i need to change it make it work for V6. 
    z=[]
	for j=1:size(y,1)
		count=0
		for mm=1:length(y[1,2])
			# if abs(x[1,2][mm]-y[j,2][mm])==2   # version 5 is one row two column. v6 is 2 row.
			if abs(x[2][mm]-y[j,2][mm])==2   # version 5 is one row two column. v6 is 2 row.
			count=count+1   # if count>30 next 
			end
		end
		t="$(x[1]) $(y[j,1]) $(count)\n"
	    z=vcat(z,t)     # append!(a,b) you can use append , push for value. 
	end
	return z
end
now()
yy=@parallel vcat for i=1:size(data2,1) 
	PT_count(data2[i,:],data1)
end
now()
of=open("count_of.txt", "w")
write(of, yy)  # $(typeof(m1))
close(of)
# awk '{if($1 != $2) print $0}' count_of.txt > tt    remove the same birds if we are using the same amount of birds for the parents test.

# -----------------method 1 for second part 
gID= Dict()
f = open("parent_gender");   # parentage gender.
for ln in eachline(f)
           blupID,sex = split(ln, r" ")
           gID[blupID] = chomp(sex)
end
close(f)

gIDFlag= []
f = open("parentGeno");    # genotype data.
for ln in eachline(f)
           ID1,geno = split(ln, r" ")
           gIDFlag=push!(gIDFlag,string(ID1))
end
close(f)

# awk '{print $1}' geno4PT.id_genofile | grep -wFf out_blupid.txt> ped_database.txt
# awk '{print $1,0,0}' childGeno > ped_database.txt
ped_ori=readdlm("ped_database.txt",' ',header=false)    # database pedigree read it to all string. you can creat your own pedigree. 

# they shold be same character. maybe not. 
genoPed=ped_ori[:,1:6]   # the fouth column is count the parents , can't be greater than 2.
genoPed[:,2]=0;genoPed[:,3]=0;genoPed[:,4]=0;genoPed[:,5]=0;genoPed[:,6]=0
# create a dict for pedigree. 
genoPedIndex=Dict()
for i=1:size(genoPed,1)
	genoPedIndex[genoPed[i,1]]=i
end

# create genoPed only using the count information.
f = open("count_of.txt");
for ln in eachline(f)
           ID1,ID2,errorC = split(ln, " ")
		   if parse(Int64,chomp(errorC)) > 19  # 550   # change here. 28
		      continue 
		   end
		   if typeof(ID2) != String
			 ID2=string(ID2)
		   end
           if typeof(ID1) != String
			 ID1=string(ID1)
		   end
		   if gID[ID2] == "M" #"1" "M"
		   genoPed[genoPedIndex[ID1],2]=ID2
		   genoPed[genoPedIndex[ID1],4]=genoPed[genoPedIndex[ID1],4]+1   # sire count
		   genoPed[genoPedIndex[ID1],6]=genoPed[genoPedIndex[ID1],6]+1   # both parents count. 
		   # write(of, "$(ID1) $(ID2) $(genderID[ID2]) $(errorC)\n")
		   end
		   if gID[ID2] == "F" # "2" "F"
		   genoPed[genoPedIndex[ID1],3]=ID2
		   genoPed[genoPedIndex[ID1],5]=genoPed[genoPedIndex[ID1],5]+1  # dam count            # dam can't great than 1
		   genoPed[genoPedIndex[ID1],6]=genoPed[genoPedIndex[ID1],6]+1  # both parenat count.   Both parents can't great than 2
		   # write(of, "$(ID1) $(ID2) $(genderID[ID2]) $(errorC)\n")
		   end
end
close(f)

# write the geno ped file. 
of11=open("geno_ped.txt", "w")
for i=1:size(genoPed,1)
write(of11, "$(genoPed[i,1]) $(genoPed[i,2]) $(genoPed[i,3]) $(genoPed[i,4]) $(genoPed[i,5]) $(genoPed[i,6])\n")  # $(typeof(m1))
# write(of11, "$(genoPed[i,1]),$(genoPed[i,2]),$(genoPed[i,3]),$(genoPed[i,4])\n")  # $(typeof(m1))
end
close(of11)


# stop here. 
grep '\.1' parentGeno


# -------------------- stop here ------------- if comparion was not needed.
# check this must. 
sort -k4n geno_ped.txt | tail -30    # check the total number of parents is here. if greater than 1, stop.
sort -k5n geno_ped.txt | tail -30     # check the total number of parents is here. if greater than 1, stop.
sort -k6n geno_ped.txt | tail -30     # check the total number of parents is here. if greater than 2, stop.

# grep P2_D10_CRYOLN_5M2F28 count_of.txt | sort -k3n | head

# parents has the same genotype. 
#sort count_of.txt | uniq > tt
#wc -l count_of.txt tt
perl /usr/local/pda/bin/getPedigree.pl BH 58 231 251 BH58_out_blupid.txt
perl /usr/local/pda/bin/getPedigree.pl GM 12 236 255 GM12_out_blupid.txt



grep F3_Hyline count_of.txt | sort -k3n | head
grep F5_Hyline count_of.txt | sort -k3n | head
grep B2_CBL12M-CBL35F count_of.txt | sort -k3n | head

# no low call samples. 
# check the number of parent found. if found more than 3. stop.



12AUM5585 6bqe9847 6BQG4874 1 2 3
12AUM5764 6bqe9861 6BQG4726 1 2 3
12AUM6945 6bqe9847 6BQG4874 2 1 3
12BQN1256 6bqe9847 6BQG4874 1 2 3
12BQO0719 6BQE9767 6BQG4962 1 2 3
12BQO0970 5bnz0234 5BNY9789 2 1 3
13BQN3121 5bnz0252 5BNY9781 2 1 3
13BQN5616 7BQG9119 7BQG9239 1 2 3
13BQN6627 6bqe9786 6BQG4974 2 1 3
16BYI7953 9AUK7758 10AUK8084 2 1 3
12BQN2866 6bqe9861 7BQG9382 2 3 5
13BQN3984 10AUK9072 11AUM3544 1 4 5
13BQN6678 11AUM4440 8BQG8929 11 6 17
16BYI8288 8BQG8403 8BQG8998 7 12 19
12BQN1956 6bqe9786 8BQG8789 41 32 73
12BQO0593 9AUK7844 9AUK7796 313 375 688
12AUM6423 9AUK7896 9AUK7870 406 541 947
12BQN2141 9AUK7887 9AUK7878 473 653 1126
12BQN1208 9AUK7820 9AUK7879 632 755 1387
12AUM5213 9AUK7839 9AUK7856 622 768 1390
12BQO0995 9AUK7896 9AUK7895 678 804 1482
12AUM5270 9AUK7885 9AUK7890 952 1138 2090
12BQO0424 9AUK7896 9AUK7890 951 1140 2091
12BQN1080 9AUK7896 9AUK7895 2714 3379 6093
13BQN3861 9AUK7896 9AUK7895 2890 3548 6438
12BQN2720 9AUK7896 9AUK7891 3179 3811 6990
13BQN5085 9AUK7896 9AUK7890 3408 4182 7590
12BQN1184 9AUK7894 9AUK7891 3561 4321 7882
12AUM6818 9AUK7896 9AUK7895 3532 4363 7895
12BQO0774 9AUK7896 9AUK7895 4044 4948 8992
12BQO0797 9AUK7896 9AUK7895 4044 4948 8992
16BYI3788 9AUK7896 9AUK7895 4044 4948 8992



Name    Call_rate
13BQN5490       0.954971857410882
13BQN3237       0.979362101313321
13BQN5259       0.943714821763602
13BQN3241       0.962476547842402
13BQN6225       0.945590994371482
13BQN5318       0.947467166979362
13BQN4369       0.966228893058161
13BQN6532       0.962476547842402
13BQN6269       0.778611632270169

==> AgriPlex_ID_results_part2.txt <==
Name    Call_rate
16BYI7496       0.874296435272045
16BYI7016       0.874296435272045
16BYI2545       0.887429643527205
16BYI5286       0.900562851782364
16BYI7292       0.866791744840525
16BYI2674       0.861163227016886
16BYI7677       0.887429643527205
16BYI2537       0.911819887429643
16BYI4647       0.857410881801126
chenju@cvirdlux05> grep -wFf tt AgriPlex_ID_results_part1.txt
13BQN5085       0.0544090056285178
13BQN3984       0.236397748592871
13BQN6678       0.270168855534709
13BQN3861       0.0825515947467167
12BQN1956       0.212007504690432
12BQN2866       0.887429643527205
12BQN2141       0.433395872420263
12BQN1184       0.0318949343339587
12AUM6423       0.714821763602251
12BQN1208       0.24015009380863
12AUM5213       0.116322701688555
12BQN2720       0.0525328330206379
12AUM5270       0.202626641651032
12BQO0424       0.326454033771107
12BQO0593       0.25703564727955
12BQO0797       NA
12BQO0774       NA
12AUM6818       0.0637898686679174
12BQO0995       0.236397748592871
12BQN1080       0.0694183864915572
chenju@cvirdlux05> grep -wFf tt AgriPlex_ID_results_part2.txt
16BYI8288       0.358348968105066
16BYI3788       NA
chenju@cvirdlux05> vi tt2
chenju@cvirdlux05> vi tt2
chenju@cvirdlux05> vi tt2
chenju@cvirdlux05> grep -wFf tt2 AgriPlex_ID_results_part2.txt
16BYI7953       0.893058161350844
chenju@cvirdlux05> grep -wFf tt2 AgriPlex_ID_results_part1.txt
13BQN5616       0.803001876172608
13BQN3121       0.919324577861163
13BQN6627       0.801125703564728
12BQO0719       0.575984990619137
12AUM5764       0.949343339587242
12AUM5585       0.943714821763602
12BQO0970       0.626641651031895
12BQN1256       0.956848030018762
12AUM6945       0.960600375234522







# fixe the -9999 parents. 
for i=1:size(genoPed,1)
	if(string(genoPed[i,2])!="0" && string(genoPed[i,3])!="0")
		continue
	end
	for j=1:size(ped_ori,1)
 	    if(string(genoPed[i,1])==string(ped_ori[j,1]))
		   println("$(genoPed[i,1])")
		   if(string(genoPed[i,2])=="0" && !(string(ped_ori[j,2]) in gIDFlag))
			  println("sire $(ped_ori[j,2]) do not have genotype!")  # write to delete list.
			  genoPed[i,2]=ped_ori[j,2]
		   end
		   if(string(genoPed[i,3])=="0" && !(string(ped_ori[j,3]) in gIDFlag))
			  println("dam $(ped_ori[j,3]) do not have genotype!")  # write to delete list.
			  genoPed[i,3]=ped_ori[j,3]
		   end
	   end
    end
end

# create the final genoped if some sires or dams is missing
# and output the both not match and sire not match.
of1=open("Delete_both.txt", "w")
of2=open("Corr_ped.txt", "w")
write(of1, "$("Blup ID"),$("Orig Sire"),$("Orig Dam"),$("Geno Sire"),$("Geno Dam"),$("Difference")\n")
write(of2, "$("Blup ID"),$("Orig Sire"),$("Orig Dam"),$("Geno Sire"),$("Geno Dam"),$("Difference")\n")
for i=1:size(genoPed,1)
	for j=1:size(ped_ori,1)
	   if(string(genoPed[i,1])==string(ped_ori[j,1]))
		   # println("$(genoPed[i,1])")
		   # if(string(genoPed[i,1])=="222AJL5469")    #  debug
		      # println("$(genoPed[i,1:3]) $(ped_ori[j,1:3])")
		   # end
		   if(string(genoPed[i,2])!=string(ped_ori[j,2]) && string(genoPed[i,3])!=string(ped_ori[j,3]))
			  println("Both parents do not match for $(genoPed[i,1])")  # write to delete list.
			  write(of1, "$(genoPed[i,1]),$(ped_ori[j,2]),$(ped_ori[j,3]),$(genoPed[i,2]),$(genoPed[i,3]),$("both")\n")
		   end
		   if(string(genoPed[i,2])!=string(ped_ori[j,2]) && string(genoPed[i,3])==string(ped_ori[j,3]))
		      println("Sire do not match for $(genoPed[i,1])")   # write to the sire list.
			  write(of2, "$(genoPed[i,1]),$(ped_ori[j,2]),$(ped_ori[j,3]),$(genoPed[i,2]),$(genoPed[i,3]),$("sire")\n")
		   end
		   if(string(genoPed[i,2])==string(ped_ori[j,2]) && string(genoPed[i,3])!=string(ped_ori[j,3]))
		      println("Dam do not match for $(genoPed[i,1])")    # write to the dam list. 
			  write(of2, "$(genoPed[i,1]),$(ped_ori[j,2]),$(ped_ori[j,3]),$(genoPed[i,2]),$(genoPed[i,3]),$("dam")\n")
		   end
	   end
    end
end
close(of1)
close(of2)
now()
# --------------- finished here. very good.
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
julia -p 30     #  30 cpus only 6 mins, 15 cpus only 15 mins.   # very good. 
# duplicate genotype checking for two files. 
# ---------------method 1.  parallele------------------------------------ 
# 
data1=readdlm("parentGeno",' ',header=false)   # k2 
data2=readdlm("childGeno",' ',header=false)    # This is for the curent MG.  k1

# data1=readdlm("no228",' ',header=false)   # k2 
# data2=readdlm("232",' ',header=false)    # This is for the curent MG.  k1

# addprocs(10)
@everywhere function PT_count(x,y)     # it is work in V5, i need to change it make it work for V6. 
    z=[]
	for j=1:size(y,1)
		count=0
		for mm=1:length(y[1,2])
			if (x[2][mm] != y[j,2][mm])   # version 5 is one row two column. v6 is 2 row.
			count=count+1   # if count>30 next 
			end
		end
		t="$(x[1]) $(y[j,1]) $(count)\n"
	    z=vcat(z,t)     # append!(a,b) you can use append , push for value. 
	end
	return z
end
now()
yy=@parallel vcat for i=1:size(data2,1) 
	PT_count(data2[i,:],data1)
end
now()
of=open("count_of_dup.txt", "w")
write(of, yy)  # $(typeof(m1))
close(of)




# pedigre check

cd /work/31/Jun/agriPlex/12232022

# awk '{print $4$5,$6$7,$8$9,$10}' PEDI_DV.txt > pedi.txt

# # checking the missing parents.
# awk '{print $6,$8}' pedigree_MG15.txt | sort | uniq -c
# # gender check.
# awk '{print $10}' pedigree_MG15.txt | sort | uniq -c

# parent gender. 
awk '{print $4$5,$6$7,$8$9,$10}' PEDI_DV.txt > pedi_parents.txt


sed -i 's/TS29_//g' geno_ped.txt
sed -i 's/BH58_//g' geno_ped.txt



library(dplyr)
# setwd("/work/31/Jun/agriPlex/Pilot3/60Kv06/PT/")
setwd("/work/31/Jun/agriPlex/12232022")
Map1=read.table(file="geno_ped.txt",header=FALSE,sep=" ",stringsAsFactors=FALSE)
Map2=read.table(file="pedi_parents.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)

colnames(Map1)=c("Animal","sire1","dam1","count_sires","count_dams","count_parents")
colnames(Map2)=c("Animal","sire2","dam2","MG","sex","hatch","source")
# colnames(Map2)=c("chr2","pos2","pos2_1","Name")

# Map2$chr2=gsub("chr","",Map2$chr2)
Map1$Animal=toupper(Map1$Animal)
Map1$sire1=toupper(Map1$sire1)
Map1$dam1=toupper(Map1$dam1)


result = inner_join(Map1,Map2,by="Animal")
dim(Map1);dim(Map2); dim(result)



any(result$sire1==result$sire2)
any(result$dam1==result$dam2)


head(result)
# all(result$chr2==result$chr1)
# all(result$pos2==result$pos1)

# Compare 
28 2 5

dim(result[which(result$sire1!=result$sire2),])
dim(result[which(result$dam1!=result$dam2),])

dim(result[which(result$dam1!=result$dam2 & result$sire1!=result$sire2),])
dim(result[which(result$dam1==result$dam2 & result$sire1!=result$sire2),])
dim(result[which(result$dam1!=result$dam2 & result$sire1==result$sire2),])

result[which(result$sire1!=result$sire2),]
result[which(result$dam1!=result$dam2),]
result[which(result$dam1!=result$dam2 & result$sire1!=result$sire2),]
result[which(result$dam1==result$dam2 & result$sire1!=result$sire2),]
result[which(result$dam1!=result$dam2 & result$sire1==result$sire2),]
? 6bqe9786  ? geneseek problem ??
dam call rate is not good. 
R> result[which(result$sire1!=result$sire2),]
       Animal    sire1      dam1     sire2      dam2 gender
131 15BYJ5986 9AUK7758 10AUK8084 10AUK8126 10AUK8084      2   ??? conflict count is 5 
138 15BYJ6011 9AUK7758 10AUK8084 10AUK8126 10AUK8084      1   ??? conflict count is 5
194 15BYJ7207 9AUK7896  9AUK7895 10AUK8348 10AUM0104      2   ??? 0 call sample
261 15BYJ9270        0  7BQG9336        00  7BQG9336      2     OK
302 15BYJ9865 9AUK7896  9AUK7895  9AUK7737  9AUK7702      2    ??? 0 call samples
R> result[which(result$dam1!=result$dam2),]
       Animal     sire1      dam1     sire2      dam2 gender
194 15BYJ7207  9AUK7896  9AUK7895 10AUK8348 10AUM0104      2  ??? 0 call samples
196 15BYJ7257 10AUK7920 10AUK9905 10AUK7920 10AUK8253      2    # why is , here is right.    8253 has big conflict.    different ID from 3K. and 26.7% and also found the different parents.
302 15BYJ9865  9AUK7896  9AUK7895  9AUK7737  9AUK7702      2    ??? 0 call samples.


15BYJ7257 is mixed samples. mother is not right. 
15BYJ9865 15BYJ7207 is 0 call sample.



grep 13BQN7408 count_of.txt | sort -nk3 | head
grep 13BQN7408 count_of.txt | sort -hk3 | head






grep 16BYI7953 count_of.txt
grep 16BYI7953 count_of.txt | sort -hk3 | head
grep 13BQN3556 count_of.txt | sort -hk3 | head   @ QC for MAP


chenju@cvirdlux04> grep 16BYI7953 count_of.txt | grep 10AUK8253
15BYJ7257 10AUK8253 54

# this bird is alos problem even it is right. 
15BYJ9806 9AUK7806 0
15BYJ9806 9AUK7882 0
15BYJ9806 9AUK7758 3



chenju@cvirdlux05> grep -wFf tt_PT AgriPlex_ID_results_part1.txt
13BQN6269       0.778611632270169
13BQN5085       0.0544090056285178
13BQN3984       0.236397748592871
13BQN6678       0.270168855534709
13BQN3861       0.0825515947467167
12BQN2141       0.433395872420263
12BQN1184       0.0318949343339587
12AUM6423       0.714821763602251
12BQN1208       0.24015009380863
12BQO0579       0.977485928705441
12AUM5213       0.116322701688555
12BQN2720       0.0525328330206379
12BQO0066       0.960600375234522
12AUM6470       0.924953095684803
12AUM5270       0.202626641651032
12BQO0424       0.326454033771107
12BQO0593       0.25703564727955
12BQO0797       NA
12BQO0774       NA
12AUM6818       0.0637898686679174
12AUM6256       0.962476547842402
12BQO0995       0.236397748592871
12BQN1080       0.0694183864915572
chenju@cvirdlux05> grep -wFf tt_PT AgriPlex_ID_results_part2.txt
16BYI3788       NA
14BQN8558       0.923076923076923
14BYJ0665       0.928705440900563
14BQN8162       0.932457786116323
14BQN8125       0.896810506566604



library(dplyr)

setwd("/work/31/Jun/rapidgenome/3152023/PT")
Map1=read.table(file="geno_Gender_part.csv",header=FALSE,sep=",",stringsAsFactors=FALSE)
Map2=read.table(file="pedi_parents.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE,skip=0)


# colnames(Map1)=c("Animal","callRate","Sex_estimate")
# colnames(Map2)=c("Animal","sire2","dam2","gender")
colnames(Map1)=c("Animal","QC1","QC2","Sex_estimate")
colnames(Map2)=c("Animal","sire2","dam2","MG","sex","hatch","source")

# Map2$chr2=gsub("chr","",Map2$chr2)
# Map1$Animal=toupper(Map1$Animal)


result = inner_join(Map1,Map2,by="Animal")
dim(Map1);dim(Map2); dim(result)


table(result$Sex_estimate,result$sex)


result[which(result$Sex_estimate == "F" & result$gender !=2),]
result[which(result$Sex_estimate == "M" & result$gender !=1),]



