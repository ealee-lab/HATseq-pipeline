chimera_list = []
nonchimera_list = []
misalignment_list = []

with open(snakemake.input.breakpoint_fa) as fa:
	peak = fa.readline().strip("\n").strip(">")
	seq = fa.readline().strip("\n").upper()
	
	while(peak):
		if "minus" in peak:
			chimera = "A"
		else:
			chimera = "T"
		
		overall_max = 0
		i = 0	
		
		while i < len(seq):
			if seq[i] == chimera:
				cur_max = 1
				i = i+1

				if i < len(seq):
					while seq[i] == chimera:
						cur_max = cur_max + 1
						i = i + 1
						if i == len(seq):
							break

				if cur_max > overall_max:
					overall_max = cur_max
			i = i+1

		if overall_max > 12 or (seq.count(chimera) / len(seq)) >= 0.8:
			chimera_list.append(peak)
		else:
			nonchimera_list.append(peak)
		peak = fa.readline().strip("\n").strip(">")
		seq = fa.readline().strip("\n").upper()

with open(snakemake.input.peak_seq) as fa:
	peak = fa.readline().strip("\n").strip(">")
	seq = fa.readline().strip("\n").upper()
	
	while(peak):
		if "minus" in peak:
			chimera = "A"
		else:
			chimera = "T"

		if (seq.count(chimera)/ len(seq) > 0.5): 
			misalignment_list.append(peak)
		peak = fa.readline().strip("\n").strip(">")
		seq = fa.readline().strip("\n").upper()
		
with open(snakemake.output.misalignment,'w') as mis:
	mis.write("\n".join(misalignment_list))		
with open(snakemake.output.chimera,'w') as chim:
	chim.write("\n".join(chimera_list))
with open(snakemake.output.nonchimera,'w') as nonchim:
	nonchim.write("\n".join(nonchimera_list))
