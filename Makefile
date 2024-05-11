name := main

result:
	mpif90 heat.f95 -o $(name)
	mpiexec -np 4 ./$(name)

clear:
	rm -f *.o $(name) data*.dat
