import pandas as pd
import gurobipy as gp
import csv
from itertools import combinations
from gurobipy import GRB
# ---- READ FILES ---- #
UB = pd.read_csv('iJO1366_UB.csv',header=None)[0].tolist()

LB = pd.read_csv('iJO1366_LB.csv',header=None)[0].tolist()

rxn = pd.read_csv('iJO1366_rxn.csv',header=None)[0].tolist()

met = pd.read_csv('iJO1366_met.csv',header=None)[0].tolist()

b = pd.read_csv('iJO1366_b.csv',header=None)[0].tolist()

c = pd.read_csv('iJO1366_c.csv',header=None)[0].tolist()

S = pd.read_csv('iJO1366_S.csv',header=None).values

r = pd.read_csv('iJO1366_r.csv',header=None)[0].tolist() # if reaction i is reversible r[i] = 1

names = pd.read_csv('iJO1366_names.csv',header=None)[0].tolist()

fba = pd.read_csv('iJO1366_FBA.csv',header=None)[0].tolist()

# --- Inner Model ----
def inner(yoj):
    i = gp.Model()
    i.params.LogToConsole = 0

    vi = i.addVars(range(len(rxn)),lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='vi')

    i.setObjective(vi[bm_core], GRB.MAXIMIZE)

    i.addConstrs((gp.quicksum(S[i,j] * vi[j] for j in range(len(rxn))) == 0 for i in range(len(met))), name='IFFBA')

    i.addConstrs((LB[j]*yoj[j] <= vi[j] for j in range(len(rxn))), name='iLB')
    i.addConstrs((vi[j] <= UB[j]*yoj[j] for j in range(len(rxn))), name='iUB')

    i.optimize()

    if i.status == GRB.OPTIMAL:
        soi = i.getAttr('X',i.getVars())
        vis = soi.copy()
    elif i.status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
        vis = yoj
        vis[bm_core] = 2000
    return vis

# --- Change biological assumptions ----

# prespecified amount of glucose uptake 10 mmol/grDW*hr 'EX_glc__D_e' = -10 reaction

LB[rxn.index('EX_glc__D_e')] = -10
UB[rxn.index('EX_glc__D_e')] = -10

# Unconstrained uptake routes for inorganic phosphate, sulfate and ammonia
# 'EX_o2_e';'EX_pi_e';'EX_so4_e'; 'EX_nh4_e' index = 184 ; 199 ; 259 ; 169

LB[rxn.index('EX_o2_e')] = 0
LB[rxn.index('EX_pi_e')] = -1000
LB[rxn.index('EX_so4_e')] = -1000
LB[rxn.index('EX_nh4_e')] = -1000

#Enable secretion routes for acetate, carbon dioxide, ethanol, formate, lactate and succinate
# {'EX_ac_e';'EX_co2_e';'EX_etoh_e';'EX_for_e';'EX_lac__D_e';'EX_succ_e'} change in the upper bound
# index = 87; 2; 340; 422; 91; 261

sec_routes = ['EX_ac_e','EX_co2_e','EX_etoh_e','EX_for_e','EX_lac__D_e','EX_succ_e']

for i in sec_routes:
    UB[(rxn.index(i))] = 1000

# Constrain the phosphotransferase system - 'b' means a change in both bounds to fix a flux
# 'GLCabcpp', -1000, 'l' - 'GLCptspp', -1000, 'l' - 'GLCabcpp', 1000, 'u' - 'GLCptspp', 1000, 'u' - 'GLCt2pp', 0, 'b'
# index 1287, 1291, 1287, 1291, 1293

LB[rxn.index('GLCabcpp')] = -1000
LB[rxn.index('GLCptspp')] = -1000
UB[rxn.index('GLCabcpp')] = 1000
UB[rxn.index('GLCptspp')] = 1000
LB[rxn.index('GLCt2pp')] = 0
UB[rxn.index('GLCt2pp')] = 0

# Defining the index for the biomass 'BIOMASS_Ec_iJO1366_core_53p95M'

bm_core = rxn.index('BIOMASS_Ec_iJO1366_core_53p95M')

# Defining the index for the chemical target succinate 'EX_succ_e'

chemical = rxn.index('EX_succ_e')
# Defining the knockout set and the combinations
set_reaction = ['GLCabcpp', 'GLCptspp', 'HEX1', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD','PGK', 'PGM', 'ENO', 'PYK',
'LDH_D', 'PFL', 'ALCD2x', 'PTAr', 'ACKr','G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TKT1', 'TALA', 'TKT2', 'FUM',
'FRD2', 'SUCOAS', 'AKGDH', 'ACONTa', 'ACONTb', 'ICDHyr', 'CS', 'MDH','MDH2', 'MDH3', 'ACALD']

knockout = []

for i in set_reaction:
    knockout.append(rxn.index(i))

comb = combinations(knockout,2)

lcomb = list(comb)

s = {}
for i in range(len(lcomb)):
    a,b = lcomb[i]
    #print('Reaction Combination:',rxn[a],'-',rxn[b])
    y = [1 for i in range(len(rxn))]
    y[a]=0
    y[b]=0
    s[(a,b)] = y

r = {}
n = 0
for i in s.keys():
    n += 1
    print('Combination %d'%n,'-',i)
    print('** Inner Model Begin**')
    vi = inner(s[i])
    if vi[bm_core] >= .5*fba[bm_core] and vi[bm_core] < 2000:
        r[i] = dict({'Biomass':vi[bm_core],'Succinate':vi[chemical]})
        print('Dictionary updated')

df = pd.DataFrame(r).T
df1 = pd.DataFrame(r)
df.to_csv('iJO1366_K2_Combinations_Suc.csv')
df1.to_csv('iJO1366_K2_Combinations_Suc_1.csv')
