## Calculation of steady-state temperature profile

#Step 1: Installation & Set up
The following pipeline merges two separate code sets. The first code (gDDA) calculates the electric fields and polarizations of a nanoparticle assembly. The second code (tDDA) uses those outputs to source the temperature calculation. 
* Compline the code "Lattice_Diffusion.c" by typing "make all" inside the source_code folder. This should make the C executable.
* Compile gDDA from my gDDA repository (https://github.com/claireannewest/gDDA). 
* In order to run tDDA, you will need a look-up table for the values of the lattice green function at all points in the temperature calculation. I have included a python script that will calculate this in the folder "lattice_greenfunction". The code is very slow, so I have included some previously calculated grids. A 20x20x20 grid is in the "lattice_greenfunction" folder. This will suffice for small structures. A much larger  300x300x300 grid can be found at: https://drive.google.com/file/d/1s4z__Sj6Ze3STPznusq1IBdATIJjLmko/view?usp=sharing. 


#Step 2: Running the code
The set of scripts that will run the above codes live in the folder "pipeline". Small shapes (10 nm radius sphere) can be run locally. For larger shapes, it is adviced to use a supercomputer.
* Design your shape using the template file "shape.f90". It is currently set up to run a 10 nm radius sphere.
* Update system parameters in the file "parameters.input". (refractive indicies, thermal conductivities, etc.)
*



To see the derivations of this code, see:
https://www.overleaf.com/read/mrdzxbrwspqt
