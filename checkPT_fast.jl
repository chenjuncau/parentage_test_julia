# Jun Chen
# --------
# awk '{print $1}' geno4PT.id_genofile | grep -wFf - out_blupid.txt> ped_database.txt
# This is much more fast, detect sire side and dam side sepearately. 
# This is working. 

# mg=HV58_266
# # awk '{print $1}' geno4PT.id_genofile | grep -wFf - out_blupid.txt> ped_database.txt
# awk '{print $1}' geno4PT.id_genofile | grep -wFf - out_blupid.txt > ped_database.txt
# grep ^${mg} geno4PT >childGeno
# grep -v ^${mg} geno4PT >parentGeno   # lowercase V
# sed -i 's/HV...//g' childGeno



# awk '{if($2==1)print $1}' parent_gender | grep -wFf - geno4PT > genoSire
# awk '{if($2==2)print $1}' parent_gender | grep -wFf - geno4PT > genoDam

# sed -i 's/HV...//g' childGeno

# julia -p 30

lenCount=815  # 815 for L12. 
dataSire=readdlm("genoSire",' ',header=false)   # k2 
dataDam=readdlm("genoDam",' ',header=false)   # k2 
data2=readdlm("childGeno",' ',header=false)    # This is for the curent MG.  k1
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
@everywhere function PT_count_one(x,y)     # it is work in V5, i need to change it make it work for V6. 
	count=0
		for mm=1:length(y[2])
			if abs(x[2][mm]-y[2][mm])==2   # 
			count=count+1   # if count>30 next 
			end
		end
		# t="$(x[1]) $(y[1]) $(count)\n"
	return count
end

now()
ped_ori=readdlm("ped_database.txt",'\t',header=false)    # database pedigree read it to all string. 

pedFlag=Dict()
for i=1:size(ped_ori,1)
	pedFlag[ped_ori[i,1]]=i
end

gIDFlag= []
for i=1:size(dataSire,1)
    gIDFlag=push!(gIDFlag,string(dataSire[i,1]))
end

for i=1:size(dataDam,1)
    gIDFlag=push!(gIDFlag,string(dataDam[i,1]))
end

# childFlag=Dict()
# for i=1:size(data2,1)
	# childFlag[data2[i,1]]=i
# end

# sire part
# -------------------------------------------------------------
# -------------------------------------------------------------
parentFlag=Dict()
for i=1:size(dataSire,1)
	parentFlag[dataSire[i,1]]=i
end

now()     # This one do not need parallele.

xxSire=[]
# size(dataSire,1)
childIndex=[]
# of=open("QchildGeno.txt", "w")
for i=1:size(data2,1)    # loop 1:0 not running. that is ok.
	# println("$(i)\n")     # something will happan if parents do not have genotype. add a if else 
	# if ( !(string(ped_ori[pedFlag[data2[i,1]],2]) in gIDFlag) || !(string(ped_ori[pedFlag[data2[i,1]],3]) in gIDFlag))
	if ( !(string(ped_ori[pedFlag[data2[i,1]],2]) in gIDFlag) )
		# write(of, "$(data2[i,1]) $(data2[i,2])\n")
         childIndex=push!(childIndex,i)		
    elseif (PT_count_one(data2[i,:],dataSire[parentFlag[ped_ori[pedFlag[data2[i,1]],2]],:]) <= lenCount)
         countSire=PT_count_one(data2[i,:],dataSire[parentFlag[ped_ori[pedFlag[data2[i,1]],2]],:])   # could be delete
         ttt1="$(data2[i,1]) $(ped_ori[pedFlag[data2[i,1]],2]) $(countSire)\n"  # could be delete
	     xxSire=vcat(xxSire,ttt1)     # could be delete
	     continue
	else
		 # write(of, "$(data2[i,1]) $(data2[i,2])\n") 
		 childIndex=push!(childIndex,i)
	end
end
		# data2 = data2[setdiff(1:end, i), :]   # parentFlag has genotype parents. 		# data22=vcat(data22,data2[i,:])
now()
# close(of)

# data22=readdlm("QchildGeno.txt",' ',header=false) 

# This only for mass mating data. mass mating data. data hsa a lot of error.use for everyting. 
yySire=[]
if(size(childIndex,1)==0)
   # error("There is no error, everyting is clean")
   println("There is no sire error, everyting is clean\n")
else
    data22=data2[childIndex,:]   # need to debug here. if childIndex==0
	yySire=@parallel vcat for i=1:size(data22,1) 
		PT_count(data22[i,:],dataSire)
	end
end

now()


# -------------------------------------------------------------
# -------------------------------------------------------------



# Dam part
# -------------------------------------------------------------
# -------------------------------------------------------------
parentFlag=Dict()
for i=1:size(dataDam,1)
	parentFlag[dataDam[i,1]]=i
end

now()     # This one do not need parallele.

xxDam=[]
# size(dataDam,1)
# of=open("QchildGeno.txt", "w")
childIndex=[]
for i=1:size(data2,1)    # loop 1:0 not running. that is ok.
	# println("$(i)\n")     # something will happan if parents do not have genotype. add a if else 
	# if ( !(string(ped_ori[pedFlag[data2[i,1]],2]) in gIDFlag) || !(string(ped_ori[pedFlag[data2[i,1]],3]) in gIDFlag))
	if ( !(string(ped_ori[pedFlag[data2[i,1]],3]) in gIDFlag))
		# write(of, "$(data2[i,1]) $(data2[i,2])\n") 
		childIndex=push!(childIndex,i)  # if no genotype, need to check with all others. 
    elseif (PT_count_one(data2[i,:],dataDam[parentFlag[ped_ori[pedFlag[data2[i,1]],3]],:]) <= lenCount)
         countDam=PT_count_one(data2[i,:],dataDam[parentFlag[ped_ori[pedFlag[data2[i,1]],3]],:])   # could be delete
         ttt2="$(data2[i,1]) $(ped_ori[pedFlag[data2[i,1]],3]) $(countDam)\n"   # could be delete
	     xxDam=vcat(xxDam,ttt2)     # could be delete
	     continue
	else
		 # write(of, "$(data2[i,1]) $(data2[i,2])\n") 
		 childIndex=push!(childIndex,i)
	end
end
		# data2 = data2[setdiff(1:end, i), :]   # parentFlag has genotype parents. 		# data22=vcat(data22,data2[i,:])
now()
# close(of)

# data22=readdlm("QchildGeno.txt",' ',header=false) 
# This only for mass mating data. mass mating data. data hsa a lot of error.use for everyting. 
# if(size(data22,1)==0)
yyDam=[]
if(size(childIndex,1)==0)
   # error("There is no error, everyting is clean")
   println("There is no dam error, everyting is clean\n")   
else
    data22=data2[childIndex,:]
	yyDam=@parallel vcat for i=1:size(data22,1) 
		PT_count(data22[i,:],dataDam)
	end
end


now()
# -------------------------------------------------------------
# -------------------------------------------------------------



of=open("conflict_count", "w")
write(of, xxSire)  # $(typeof(m1))  no problem parents.
write(of, yySire)  # $(typeof(m1))
write(of, xxDam)  # $(typeof(m1))  no problem parents.
write(of, yyDam)  # $(typeof(m1))
close(of)

# stop here. 
# stop here. 
# stop here. 
# stop here. 
# stop here. 


# # 11:48
# dir_scpt="/biotech/Jun/programs/script_Jun/"
# dir_work="/biotech/genomic_selection/HV-58/run/mg266/PT_2023_09_12_13_15"
# perl ${dir_scpt}/create_G_Pedigree_thr20.pl \
	# ${dir_work}/conflict_count \
	# ${dir_work}/geno4PT.id_genofile \
	# ${dir_work}/parent_gender \
	# ${dir_work}/ped_genomic \
	# > log_create_G_Pedigree_thr18 
	


# hdr = vec(readdlm("header.txt",'\t'))
# dat2 = dat1[:, hdr]    # hdr is vector
# dat2 = select(df1, hdr)

	
