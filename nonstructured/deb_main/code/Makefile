compiler = ifort -check -traceback
#compiler = ifort -fast -parallel
directory=../
core = prmtr.o grid.o condition.o scheme.o OUTPUT.o euler_exp.o
all: a.out

prmtr.o:prmtr.f90
	$(compiler) -c $<
grid.o: grid.f90 prmtr.f90
	$(compiler) -c $<
condition.o: condition.f90 prmtr.f90
	$(compiler) -c $<
scheme.o: scheme.f90 prmtr.f90
	$(compiler) -c $<
euler_exp.o: euler_exp.f90 prmtr.f90
	$(compiler) -c $<
OUTPUT.o: OUTPUT.f90 prmtr.f90
	$(compiler) -c $<

a.out: main.f90 $(core)
	$(compiler) $^

clean:remove
	rm -rf *.mod a.out *.o *.vtk result/*.vtk *.bin *.d
remove:
mesh_clean:
	rm -rf data/*.msh
