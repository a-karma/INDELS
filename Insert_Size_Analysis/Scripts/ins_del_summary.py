import csv
import sys
import statistics as stat

vcf_file_input=sys.argv[1]
indiv_file_input=sys.argv[2]
min_indel_size=int(sys.argv[3])

with open(indiv_file_input,'r') as idv:
	idv_content=idv.readlines()
indiv_list=[el.strip() for el in idv_content]
insertion_dict=dict.fromkeys(indiv_list,[])
deletion_dict=dict.fromkeys(indiv_list,[])
event_dict=dict.fromkeys(indiv_list,[0,0,0,0,0])

with open(vcf_file_input,'r') as f:
	for line in csv.reader(f, dialect="excel", delimiter='\t'):
		if line[0][0]!='#':
			L_list=[]
			inds_info=line[9:]
			l_base=len(line[3])
			L_list.append(l_base)
			alternatives=line[4].split(",")
			for element in alternatives:
				l_element=len(element)
				L_list.append(l_element)
			for i in range(len(inds_info)):
				S_ID=indiv_list[i]
				full_geno_info=inds_info[i].split(":")
				haplo=full_geno_info[0].replace("/",";")
				haplo_2=haplo.replace("|",";")
				haplo_list=haplo_2.split(";")
				event_line=[0,0,0,0,0]
				if haplo_list[0]==".":
					event_line[0]+=1
				else:
					haplo1_len=L_list[int(haplo_list[0])]-l_base
					if haplo1_len>min_indel_size:
						insertion_dict[S_ID]=insertion_dict.get(S_ID)+[haplo1_len]
						event_line[3]+=1
					elif haplo1_len<-min_indel_size:
						deletion_dict[S_ID]=deletion_dict.get(S_ID)+[haplo1_len]
						event_line[4]+=1
					elif haplo1_len==0:
						event_line[1]+=1
					else:
						event_line[2]+=1
				if (len(haplo_list)==1 or haplo_list[1]=="."):
					event_line[0]+=1
				else:
					haplo2_len=L_list[int(haplo_list[1])]-l_base
					if haplo2_len>min_indel_size:
						insertion_dict[S_ID]=insertion_dict.get(S_ID)+[haplo2_len]
						event_line[3]+=1
					elif haplo2_len<-min_indel_size:
						deletion_dict[S_ID]=deletion_dict.get(S_ID)+[haplo2_len]
						event_line[4]+=1
					elif haplo2_len==0:
						event_line[1]+=1
					else:
						event_line[2]+=1
				event_dict[S_ID]=[sum(x) for x in zip(event_dict.get(S_ID),event_line)]

for individual in indiv_list:
	if len(insertion_dict.get(individual))>=1:
		mean_ins=round(stat.mean(insertion_dict.get(individual)),2)
	else:
		mean_ins=0.0
	if len(insertion_dict.get(individual))>=2:
		sd_ins=round(stat.stdev(insertion_dict.get(individual)),2)
	else:
		sd_ins=0.0
	if len(deletion_dict.get(individual))>=1:
		mean_del=round(stat.mean(deletion_dict.get(individual)),2)
	else:
		mean_del=0.0
	if len(deletion_dict.get(individual))>=2:
		sd_del=round(stat.stdev(deletion_dict.get(individual)),2)
	else:
		sd_del=0.0
	net_len=sum(insertion_dict.get(individual)+deletion_dict.get(individual))
	event_dict[individual]=[str(value) for value in event_dict.get(individual)+[mean_ins,sd_ins,mean_del,sd_del,net_len]]

outname="Output/Indels_summary_minSize_"+str(min_indel_size)+".tsv"

with open(outname, 'w') as osf:
	for val in event_dict.keys():
		summary="\t".join(event_dict.get(val))
		outline=val+"\t"+summary+"\n"
		osf.write(outline)
