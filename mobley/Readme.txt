Directions for running the Mobley_Surface_Area, Param_Water_With_Ions, Param_With_Outliers, and Random_Test_Set branches.

Mobley_Surface_Area:

This branch tests the models efficacy in predicting solvation free energies for a consistent training set across all 9 selected solvents.  In order to run an optimization for a particular solvent, open ParamSolvents.py.  

1) Modify the solvents you want to optimize in the list under <if RunType == 'Param'>.  For example, if you want to parameterize for a solvent not in the list, simply add it to the list.  

2) You can select specific entries from the list by slicing the Solvent_List in line 16.  e.g. if you want the last solvent only set Solvent_List = Solvent_List[-1].

3) In the terminal type python ParamSolvents.py.

4) To get results, run the loadResults.m file in Matlab.  This will plot all of the calculated vs. predicted free energy plots, histograms of errors, and will calculate all of the rms errors and other metrics we want.

Param_Water_With_Ions:

This branch works the same as the Mobley_Surface_Area branch except, there are three ions in the training set for water.

Param_With_Outliers:

This branch works the same as the Mobley_Surface_Area branch except, for each solvent's training set, we have added the two solutes that yielded the largest errors between calculated free energies and reference free energies.

Random_Test_Set:

This branch works the same as the Mobley_Surface_Area branch except, a random training set of 12 solutes is selected to parameterize each solvent model.  These solutes are taken from the list of 27 solutes in common to all 9 solvents.
