In order to compile a code named: mycode.f08 with a module named: mymodule.f08 and compiled as: myprogram.out run the following terminal commands.

gfortran -c mymodule.f08
gfortran mycode.f08 mymodule.o -o myprogram

The code can be compiled with more modules like so

gfortran mycode.f08 mymodule1.o mymodule2.o [...] mymoduleN.o

In order to check for uninitialised and unused variables run

gfortran -g -fcheck=all -Wall mycode.f08

sc = simple avoided crossing
dc = double avoided crossing
ec = extended coupling
sb = spin-boson
pes = potential energy surface evaluation
traj = trajectories (this was used to write the trajectories used in the thesis as examples)
