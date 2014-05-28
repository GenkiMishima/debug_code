compiler=ifort -traceback -check
#compiler=ifort -xHOST -O3 -static -parallel -openmp
compiler1=ifort -xHOST -O3 -static

core = n_grid.o prmtr.o grid.o sub_chem.o func_chem.o mod_chem.o variable.o condition.o OUTPUT.o muscl.o lusgs.o euler_exp.o scheme.o viscous.o

all: a.out

n_grid.o: n_grid.f90
	$(compiler) -c $<
prmtr.o: prmtr.f90
	$(compiler) -c $<
mod_chem.o: mod_chem.f90 n_grid.o
	$(compiler) -c $<
func_chem.o: func_chem.f90 mod_chem.o
	$(compiler) -c $<
sub_chem.o: sub_chem.f90 mod_chem.o func_chem.o
	$(compiler) -c $<
muscl.o: muscl.f90
	$(compiler) -c $<
variable.o: variable.f90
	$(compiler) -c $<
grid.o: grid.f90 variable.o
	$(compiler) -c $<
condition.o: condition.f90 variable.o
	$(compiler) -c $<
lusgs.o: lusgs.f90 variable.o
	$(compiler) -c $<
euler_exp.o: euler_exp.f90 variable.o
	$(compiler) -c $<
OUTPUT.o: OUTPUT.f90 variable.o
	$(compiler1) -c $<
scheme.o: scheme.f90 variable.o
	$(compiler) -c $<
viscous.o: viscous.f90 variable.o
	$(compiler) -c $<

a.out: main.f90 $(core)
	 $(compiler) $^
#	ifort 
#	ifort -fast -parallel 
#  ifort -traceback -idb -check all
clean:remove
	rm -rf *.mod a.out *.o data/*.d *.d *.bin
remove:
tar:
	tar zcvf data.tar.gz data/*

