#compiler = ifort -check -traceback
compiler = ifort -fast -parallel
directory=../
core = prmtr.o grid.o
all: a.out

prmtr.o:prmtr.f90
	$(compiler) -c $<
grid.o: grid.f90 prmtr.f90
	$(compiler) -c $<

a.out: postprocess.f90 $(core)
	$(compiler) $^

clean:remove
	rm -rf *.mod a.out *.o *.vtk result/*.vtk *.bin *.d
remove:
mesh_clean:
	rm -rf data/*.msh
