at_list = []

with open(snakemake.input.peak_seq) as fa:
	peak = fa.readline().strip("\n").strip(">")
	seq = fa.readline().strip("\n").upper()
	
	while(peak):
		if "minus" in peak:
			chimera = "A"
		else:
			chimera = "T"

		chimera_prop = seq.count(chimera) / len(seq)
		at_list.append(f"{peak}\t{chimera_prop}")
		
		peak = fa.readline().strip("\n").strip(">")
		seq = fa.readline().strip("\n").upper()
		
with open(snakemake.output.A_T_content,'w') as at:
	at.write("\n".join(at_list))			
