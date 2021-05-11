# Libraries to import
import pandas as pd
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import csv
from itertools import combinations

# Read the data from the csv files
UB = pd.read_csv('iJO1366_UB.csv',header=None)[0].tolist()
LB = pd.read_csv('iJO1366_LB.csv',header=None)[0].tolist()
rxn = pd.read_csv('iJO1366_rxn.csv',header=None)[0].tolist()
met = pd.read_csv('iJO1366_met.csv',header=None)[0].tolist()
b = pd.read_csv('iJO1366_b.csv',header=None)[0].tolist()
c = pd.read_csv('iJO1366_c.csv',header=None)[0].tolist()
S = pd.read_csv('iJO1366_S.csv',header=None).values
r = pd.read_csv('iJO1366_r.csv',header=None)[0].tolist() # if reaction i is reversible r[i] = 1
names = pd.read_csv('iJO1366_names.csv',header=None)[0].tolist()
vjc = pd.read_csv('iJO1366_vj_COBRA_sol.csv',header=None)[0].tolist()
fba = pd.read_csv('iJO1366_FBA.csv',header=None)[0].tolist()

# Changing the biological assumptions

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

# List of set of reactions for Knockouts ie. only reactions in this set will be deleted

set_reaction = ['GLCabcpp', 'GLCptspp', 'HEX1', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD','PGK', 'PGM', 'ENO', 'PYK', 
'LDH_D', 'PFL', 'ALCD2x', 'PTAr', 'ACKr','G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TKT1', 'TALA', 'TKT2', 'FUM',
'FRD2', 'SUCOAS', 'AKGDH', 'ACONTa', 'ACONTb', 'ICDHyr', 'CS', 'MDH','MDH2', 'MDH3', 'ACALD']

knockout = []

for i in set_reaction:
    knockout.append(rxn.index(i))

# succinate production and growth rate

rxn_target = ['EX_succ_e','EX_etoh_e','EX_for_e','EX_lac__D_e','EX_ac_e']
in_target = []

for i in rxn_target:
    in_target.append(rxn.index(i))
# Defining the index for the biomass 'BIOMASS_Ec_iJO1366_core_53p95M'

bm_core = rxn.index('BIOMASS_Ec_iJO1366_core_53p95M')

# Defining the index for the chemical target succinate 'EX_succ_e'

chemical = rxn.index('EX_succ_e')

k1 = rxn.index('PFL')
k2 = rxn.index('TALA')

FBA_WT = fba[bm_core]

k = 2
reactions = [i for i in range(len(rxn))]
metabolites = [i for i in range(len(met))]

def inner(yoj):
    global vis

    i = gp.Model()
    i.params.LogToConsole = 0

    vi = i.addVars(range(len(rxn)),lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='vi')
    
    i.setObjective(2000*vi[bm_core] + vi[chemical], GRB.MAXIMIZE) 

    i.addConstrs((gp.quicksum(S[i,j] * vi[j] for j in reactions) == 0 for i in metabolites), name='IFFBA')

    i.addConstrs((LB[j]*yoj[j] <= vi[j] for j in reactions), name='iLB')
    i.addConstrs((vi[j] <= UB[j]*yoj[j] for j in reactions), name='iUB')
    
    i.optimize()

    if i.status == GRB.OPTIMAL:
        soi = i.getAttr('X',i.getVars())
        vis = soi.copy()
        
    elif i.status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
        vis = yoj
        vis[bm_core] = 2000
        
    return vis
  def lazyctr(model,where):
    #print(str(where == GRB.Callback.MIPSOL),'-','MIPSOL')
    if where == GRB.Callback.MIPSOL:
        print('** Begin Inner Model **')
        model._voj = model.cbGetSolution(model._vars) #solutions
        model._yoj = model.cbGetSolution(model._varsy)
        model._m   = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
   
        keys = model._vars.keys()  #keys
        eps = .00001
        print('** Best Obj Bnd ***')
        print(model._m)
        print('*** yj values <= (1-eps) ***')
        for i in keys:
            if model._yoj[i]<=(1-eps):    
                print('y_%d'%i,';',model._yoj[i])
            
        model._vij = inner(model._yoj)
        v_o = [i for i in range(len(rxn)) if abs(model._vij[i]) == 0]
        comb_v = list(combinations(v_o,2))
        print('   ','*** Deletion Strategy***')
        for i in reactions:
            if model._yoj[i] < .5:
                print('      ','Reaction:','->',rxn[i])
        print('   ','*** Begin Lazy Constraints ***')
        if abs(model._vij[bm_core]-model._voj[bm_core]) >= eps:
            
            print('      ','**Input Values for lazy constraints **')
            print('         ','Inner biomas:',model._vij[bm_core])
            print('         ','Outer biomas:',model._voj[bm_core])
            print('         ','Inner succinate:',model._vij[chemical])
            print('         ','Outer succinate:',model._voj[chemical])
            
            model.cbLazy(model._vij[bm_core] <= model._vars[bm_core] + 
                         model._vij[bm_core]*(gp.quicksum((1-model._varsy[j]) for j in keys if abs(model._vij[j]) >= eps)))
            for i,j in comb_v:
                model.cbLazy(model._vars[chemical] <= model._vij[chemical] + (model._m)*(model._varsy[i]+model._varsy[j]))
            
            print('   ','*** End Lazy Constraints***')
            
    elif where == GRB.Callback.MIPNODE: # try to indent it at the same level as the first "if"
        print(str(where==GRB.Callback.MIPNODE),'-','MIPNODE')
        print(str(GRB.Callback.MIPNODE_STATUS),'-','MIPNODE Status')
        print('*** Check 2 Set Solution  ***')    # try to check the status of MIPNODE
        model.cbSetSolution(model._vars, model._vij)
        model.cbSetSolution(model._varsy, model._yoj)
        print('***Node Relaxation***')
        model._m1  = model.cbGet(GRB.Callback.MIPNODE_OBJBND)
        print(model._m1)
        keys1 = model._vars.keys()
        for i in keys1:
            if abs(model.cbGetNodeRel(model._varsy)[i]) <= (1 - .00001) :
                print('y_%d:'%i,model.cbGetNodeRel(model._varsy)[i])
        
        print('*** Set Solution passed ***')
    
m = gp.Model()

v = m.addVars(range(len(rxn)), lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name= 'v')
y = m.addVars(range(len(rxn)), vtype=GRB.BINARY, name= 'y')

m.setObjective(1*v[chemical], GRB.MAXIMIZE)

m.addConstrs((gp.quicksum(S[i,j] * v[j] for j in reactions) == 0 for i in metabolites), name='OFFBA')

# Bounds and y's
m.addConstrs((y[j] == 1 for j in reactions if j not in knockout))

m.addConstrs((LB[j]*y[j] <= v[j] for j in reactions), name='LBy')
m.addConstrs((v[j] <= UB[j]*y[j] for j in reactions), name='UBy')

m.addConstr(sum(1-y[j] for j in knockout) == k, name='Knapsack') 

m.addConstr(v[bm_core] >= .5*FBA_WT, name='expected_growth') # the value for FBA_WT comes from the first FBA



m._vars = v
m._varsy = y
m.Params.lazyConstraints = 1
m.Params.LogToConsole = 0
m.Params.OutputFlag = 1
m.Params.LogFile = 'iJO1366_LOGFILE.log'
m.optimize(lazyctr)


s = m.Runtime

if m.status == GRB.OPTIMAL:
    sol = m.getAttr('X',m.getVars())
    sol1= sol.copy()
    yoj = sol1[len(rxn):]
    voj = sol1[:len(rxn)]
    
      
if m.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
    m.computeIIS()
    if m.IISMinimal:
        print('IIS is minimal\n')
    else:
        print('IIS is not minimal\n')
        print('\nThe folllowing constraint(s) cannot be satisfied:')
    for c in m.getConstrs():
        if c.IISConstr:
            print('%s'%c.constrName)

m = s/60

h = m/60

print('Run time in seconds:',' ', s)
print('Run time in minutes:',' ', m)
print('Run time in hours:',' ',h)

