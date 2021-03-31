# iJO1366-Bilevel when k = 2
This code is to estimate the knock out strategy for succinate overproduction over the iJO1366 E.Coli metabolic network. It is written in Python and uses Gurobi solver.

The metabolic network has 1805 metabolites and 2583 reactions. Thus the Stoichiometric matrix is 1805x2853, in other words the model has 2583 continuous variables and 2583 binary variables. 

# Code description

Read the csv files (provided) to read the data the code will need, pass each item as a list or a matrix, the case of 'S' (Stoichimetric Matrix) for the easy manipulation.

Manipulate the Lower and Upper bounds according to metabolic assumptions.

Define the inner model function, this takes the y values from the outer model and calculates the biomass inner value

Define the lazy constraint function, takes the outer variables and values and starts the gurobi callbacks

Define the outer model and provides an initial solution for the inner model and the callback


