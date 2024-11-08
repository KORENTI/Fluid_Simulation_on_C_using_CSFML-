# Fluid_Simulation_on_C_using_CSFML-
Implementation of fluid simulation on C programming and using CSFML liberary for 2D display 
The code is a modification of the original implementation by Mike Ash and you can check how it works on 
https://mikeash.com/pyblog/fluid-simulation-for-dummies.html and also the java script implementation on the coding train 
@https://mikeash.com/pyblog/fluid-simulation-for-dummies.html

=> you need to install C bindings for the SFML graphics display that can be found on CFML website 
   @ https://www.sfml-dev.org/download/csfml/

=> Compile using gcc and use the right flags to incorporate the libberaries of csfml 
  
 gcc Fluid01.c -o Fluid01 -lcsfml-window -lcsfml-system -lcsfml-graphics -lm 

=> run it 

./Fluid01

