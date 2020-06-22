CC=gcc
CFLAGS=-c -Wall -iquote/gscratch/chem/masiell_group/tDDA_0919_sub_plane/include -iquote/gscratch/chem/lib/include

all: TDDA

TDDA: Lattice_Diffusion.o find_neighbor.o
	$(CC) Lattice_Diffusion.o find_neighbor.o -o Lattice_Diffusion -L/gscratch/chem/lib/ -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgfortran -llapack -lblas -lm

Lattice_Diffusion.o: Lattice_Diffusion.c
	$(CC) $(CFLAGS) Lattice_Diffusion.c

find_neighbor.o: find_neighbor.c
	$(CC) $(CFLAGS) find_neighbor.c

clean:
	rm *.o Lattice_Diffusion
