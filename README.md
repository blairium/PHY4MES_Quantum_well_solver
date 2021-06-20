## PHY4MES_Quantum_well_solver

These matlab scripts make up the second computational lab for PHY4MES. These are largely the work of Alex Schenk with my own modifications.

Any errors are my own.


### Notes on completing this assignment:
The biggest challenge for this assignment, outside of understanding the code is how best to present the data. While some plotting is included 
plotting every possible plot for every possible combination of parameters would lead to an unrealsonable number of plots and make comparisons 
needlessly difficult. As such it is nesscessary to carefully consider how you display your data. Subplots, hold on and For loops are very much your friends.

### Part A:
The file Q1.ipynb is a jupyter notebook that plots the potentials for a 1D potential well. As Part A was completed entirely by me it was done is python, a superior programming language.
If you don't want to install python and jupyter notebook, consider opening a jupyter notebook on [SWAN](https://support.aarnet.edu.au/hc/en-us/sections/360000129695-CloudStor-SWAN) which is provided by [Cloudstor](cloudstor.aarnet.edu.au/) which is free for all students and researchers in Australia

### Part B
The file QW_1.m is a matlab script that iteratively solves the coupled Poisson-Schrodinger equations to simulate a single finite potential well

### Part C 
The file QW_2.m performs in the same method as QW_1.m but for two finite potential wells

### Part D
GaAs_solver.m works as above to sinulate two quantum wells in GaAs/AlGaAs