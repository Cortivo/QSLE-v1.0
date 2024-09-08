include ../../Makefile.in

# override for debugging
#COMPOPT=-g -O0 -DDEBUG

libzmat.a: atom.o frameRT.o molecule.o config.o jacobi.o matrix3D.o matrix4D.o
	ar qc libzmat.a ./*.o
	ranlib libzmat.a

atom.o: atom.cpp ./include/atom.h
	$(CPP) $(COMPOPT) $(II) -c atom.cpp $(LIBS)
frameRT.o: frameRT.cpp ./include/frameRT.h
	$(CPP) $(COMPOPT) $(II) -c frameRT.cpp $(LIBS)
molecule.o: molecule.cpp ./include/molecule.h 
	$(CPP) $(COMPOPT) $(II) -c molecule.cpp $(LIBS)
config.o: config.cpp ./include/config.h
	$(CPP) $(COMPOPT) $(II) -c config.cpp $(LIBS)
jacobi.o: jacobi.cpp ./include/jacobi.h
	$(CPP) $(COMPOPT) $(II) -c jacobi.cpp $(LIBS)
matrix3D.o: matrix3D.cpp ./include/matrix3D.h
	$(CPP) $(COMPOPT) $(II) -c matrix3D.cpp $(LIBS)
matrix4D.o: matrix4D.cpp ./include/matrix4D.h
	$(CPP) $(COMPOPT) $(II) -c matrix4D.cpp $(LIBS)


clean:
	rm -rf ./*.o ./libzmat.a
