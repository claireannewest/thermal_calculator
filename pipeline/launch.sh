#!/bin/bash
gfortran shape.f90 && ./a.out
python input_generator.py -d parameters.input shape.dat
python input_generator.py -v parameters.input shape.dat
wait
../gDDA_source_code/ddscat &> DDA.out &
# wait
# mv tdda_input_w000_ddscat.par tdda_input
# wait

# # tDDA calculation
# /gscratch/chem/masiello_group/tDDA_1123_sub_plane/Lattice_Diffusion /gscratch/chem/masiello_group/myGreen.num_300 var.par tdda_input temp.out
