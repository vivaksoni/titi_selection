import sys
import pandas as pd
import math
import os
import numpy as np
import msprime
import tskit
import libsequence
import allel
import argparse


#Example input
#python aye_aye_null_msprime.py -region 1 -seq_len 50000 -num_replicates 100 \
#-outPath "/home/vivak/"

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-region', dest = 'region', action='store', nargs = 1, type = int, help = 'Region')
parser.add_argument('-num_replicates', dest = 'num_replicates', action='store', nargs = 1, type = int, help = 'Number of replicates')
parser.add_argument('-m_map', dest = 'm_map', action='store', nargs = 1, type = str, help = 'mutation rate map')
parser.add_argument('-r_map', dest = 'r_map', action='store', nargs = 1, type = str, help = 'recombination rate map')

args = parser.parse_args()
region = args.region[0]
num_replicates = args.num_replicates[0]
m_map = args.m_map[0]
r_map = args.r_map[0]

#Demog lengths
chroms = {
1:81997209,
2:57596708,
3:59482194,
4:51452482,
5:56421937,
6:37428171,
7:57377548,
8:42072168,
9:32858669,
10:16254630,
11:13862685,
12:78475508,
13:36722319,
14:41001960,
15:32326613,
16:37251058,
17:16847824,
18:37884098,
19:19023431,
20:37470834,
21:27914525,
22:13844647
}

def titi_demog(seq_len, num_replicates, m_map, r_map):
    N_current = 12337.5
    t1 = 3160
    N1 = 1999795.5
    t2 = 131015
    N2 = 45656
    t3 = 351168
    N3 = 169078.5
    Nancestral = N3  
    
    demography = msprime.Demography()
    demography.add_population(
        name="titi",
        description="titi population",
        initial_size=N_current,
        growth_rate=0,
        default_sampling_time=0, 
        initially_active=True,
    )
    
    #Add events
    demography.add_population_parameters_change(initial_size=N1, time=t1, growth_rate=0, population="titi")
    demography.add_population_parameters_change(initial_size=N2, time=t2, growth_rate=0, population="titi")
    demography.add_population_parameters_change(initial_size=N3, time=t3, growth_rate=0, population="titi")
    demography.sort_events

    ancestry_reps = msprime.sim_ancestry(
        {"titi": 6}, 
        demography=demography, 
        sequence_length = seq_len,
        recombination_rate = r_map,
        num_replicates=num_replicates)

    for ts in ancestry_reps:
        mutated_ts = msprime.sim_mutations(ts, rate=m_map)
        yield mutated_ts


coords = [x for x in range(0, chroms[region], 1000)]
coords.append(chroms[region])
m = list(pd.read_csv(m_map, names=['m'], sep='\t')['m'])
r = list(pd.read_csv(r_map, names=['r'], sep='\t')['r'])
m_map = msprime.RateMap(position=coords, rate=m)
r_map = msprime.RateMap(position=coords, rate=r)

ts_list = []
for replicate_index, ts in enumerate(titi_demog(chroms[region], num_replicates, m_map, r_map)):
    ts_list.append(ts)

for rep, ts in enumerate(ts_list):
    lst = []
    #Loop through variants
    for var in ts.variants():
        #Obtain counts of ancestral and derived alleles
        unique, counts = np.unique(var.genotypes, return_counts=True)
        #Append list of positions and MACs to main list
        lst.append([int(var.site.position), counts.min()])
    #Convert to df
    df = pd.DataFrame(lst, columns=['position', 'x'])
    df['n'] = 12
    df['folded'] = 1
    df = df[df.x<12]

    df2 = pd.DataFrame(lst, columns=['physPos', 'x'])
    df2['n'] = 12
    df2['genPos'] = 'NA'
    df2 = df2[['physPos', 'genPos', 'x', 'n']]
    df2 = df2[df2.x<12]

    df.to_csv("chr" + str(region) + "_rep" + str(rep) + "_SF2.aff", header=True, index=False, sep='\t')
    df2.to_csv("chr" + str(region) + "_rep" + str(rep) + "_BM.aff", header=True, index=False, sep='\t')

