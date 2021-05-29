# CSECGEFoam
Continuum scale electrochemical gas evolution solver on OpenFOAM

# Overview
CSECGEFoam is Volume of Fluid (VOF) framework to simulate the evolution of continuum scale hydrogen bubble evolution during water electrolysis. The model takes into account interfacial mass transfer for bubble growth due to supersaturation, electrochemical reactions, transport of dissolved gases and charge transport. Inorder to simulate submillimeter sized bubbles Sharp Surface force model is implemented.

# Description
This version is developed for OpenFOAM 6. Please follow the following step to compile the solver:

Step 1: Go to the directory ```$ cd Libraries/SSF_transportModels/``` and run ``` $ ./Allwclean ``` followed by ``` $ ./Allwmake ```. 

Step 2: Now go the the CSECGEFoam folder ```$ cd ../../CSECGEFoam/``` and run ``` $ wclean ``` followed by ``` $ wmake ```.

This should have installed CSECGEFoam solver in ```FOAM_USER_APPBIN ```.

# Verification 
The verification of the solver is discussed in 

# Usage
The CSECGEFoam solver is released under Creative Commons Attribution-NonCommercial 4.0 International Public License (see LICENSE.txt included included in the main directory). Please feel free to extend, modify and use the solver according to your requirement (only for noncommerical usage). For reference of the solver, please cite to the following papers [1,2,3].

# Main references
1. Vachaparambil, K.J., Einarsrud, K.E. 2021. Numerical simulation of continuum scale electrochemical hydrogen bubble evolution. Applied Mathematical Modelling. [doi:10.1016/j.apm.2021.05.007(https://doi.org/10.1016/j.apm.2021.05.007)
 
2. Vachaparambil, K.J., Einarsrud, K.E. 2020. On Modelling Electrochemical Gas Evolution Using The Volume Of Fluid Method. Proceedings from the 14th International Conference on CFD in Oil & Gas, Metallurgical and Process Industries (CFD2020), SINTEF Academic Press, pp. 17-27. [https://hdl.handle.net/11250/2720838](https://hdl.handle.net/11250/2720838)

3. Vachaparambil, K.J. 2020. Interface resolved simulations of continuum scale electrochemical hydrogen evolution. Doctoral theses at NTNU, 363. [https://hdl.handle.net/11250/2723231](https://hdl.handle.net/11250/2723231)

