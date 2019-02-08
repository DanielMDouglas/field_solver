solver: solver.cpp field.cpp
	g++ solver.cpp field.cpp -o solver

clean:
	rm solver
	rm *.dat
