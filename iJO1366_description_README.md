# iJO1366-Bilevel when k = 2
This code is to estimate the knock out strategy (two reactions to be deleted from the network) for succinate overproduction over the iJO1366 E.Coli metabolic network. It is written in Python and uses Gurobi solver.

This is a bi-level problem where the outer objective is to optimize succinate overproduction (chemical of interest) and the inner objective is maximize cellular growth (biomass).

The metabolic network has 1805 metabolites and 2583 reactions. Thus the Stoichiometric matrix is 1805x2853. The outer model has 2583 continuous variables and 2583 binary variables. And, the inner model has 2583 continuous variables

Code inputs:

      UB    = Upper bounds, same lenght as rxn
      LB    = Lower bounds, same lenght as rxn
      rxn   = Reaction abreviations, a list of 2583 elements 
      met   = Metabolites abreviations, a list of 1805 elements
      S     = Stoichiometric Matrix, an array of (1805,2583)
      fba   = iJO1366 Flux Balance Analysis, reactions fluxes without the knockouts, full metabolic network

# Code description

* Read the csv files (provided) to read the data the code will need, pass each item as a list or a matrix, the case of 'S' (Stoichimetric Matrix) for the easy manipulation.

* Manipulate the Lower and Upper bounds according to metabolic assumptions.
   * This is done by accesing the indeces of the upper and lower bounds according to the reaction abreviation, example to set the glucose intake equal to 10:
   
      LB[rxn.index('EX_glc__D_e')] = -10;
      UB[rxn.index('EX_glc__D_e')] = -10
      
   * Note that the name in parentheses is the abrreviation, list.index('name') returns the index of that element in the list, since UB,LB and rxn are the same lenght the UB,LB value for that index corresponds to the reaction upper and lower bound, so it can be manipulated.
   
* Define the inner model function, this takes the y values from the outer model and calculates the biomass inner value

* Define the lazy constraint function, takes the outer variables and values and starts the gurobi callbacks

* Define the outer model and provides an initial solution for the inner model and the callback

Each section in the code is identified by '#' followed by its name
