chimera_list = []

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

		chimera_prop = seq.count(chimera) / len(seq)
		chimera_list.append(f"{peak}\t{overall_max}\t{chimera_prop}")

		peak = fa.readline().strip("\n").strip(">")
		seq = fa.readline().strip("\n").upper()

with open(snakemake.output.chimera,'w') as chim:
	chim.write("\n".join(chimera_list))
