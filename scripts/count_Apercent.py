#!/usr/bin/python 
import sys
import pandas as pd 
import numpy as np 

plus_fa = sys.argv[1]
minus_fa = sys.argv[2] 
chimera = plus_fa.replace("plus_breakpoints.fa", "chimera.txt") 
nonchimera = plus_fa.replace("plus_breakpoints.fa", "nonchimera.txt")
#plus_bed = plus_fa.strip(".fa") + ".bed" 
#minus_bed = minus_fa.strip(".fa") + ".bed"
chimera_list = []
nonchimera_list = []

with open(plus_fa) as fa:
    peak = fa.readline().strip("\n").strip(">") 
    seq = fa.readline().strip("\n")
    while(peak): 
        percent_3prime = int(seq[int(15):len(seq)-1].count("T")) / len(seq[int(15):len(seq)-1])
        if(percent_3prime < .8):
            nonchimera_list.append(peak) 
        else:
            chimera_list.append(peak)
        peak = fa.readline().strip("\n").strip(">")
        seq = fa.readline().strip("\n")

with open(plus_fa) as fa:
    peak = fa.readline().strip("\n").strip(">")
    seq = fa.readline().strip("\n")
    while(peak):
        percent_3prime = int(seq[int(15):len(seq)-1].count("T")) / len(seq[int(15):len(seq)-1])
        if(percent_3prime < .8):
            nonchimera_list.append(peak)
        else:
            chimera_list.append(peak)
        peak = fa.readline().strip("\n").strip(">")
        seq = fa.readline().strip("\n")

with open(chimera, "w") as f:
    f.write('\n'.join(chimera_list))
with open(nonchimera, "w") as f:
    f.write('\n'.join(nonchimera_list))
