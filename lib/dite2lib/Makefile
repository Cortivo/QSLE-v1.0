include ../../Makefile.in

# override for debugging
#COMPOPT=-g -O0 -DDEBUG

libdite2.a: DiTe2.o inputParser.o
	ar qc libdite2.a ./*.o
	ranlib libdite2.a

DiTe2.o: DiTe2.cpp ./include/DiTe2.h
	$(CPP) $(COMPOPT) $(II) $(LL) -c DiTe2.cpp $(LIBS) 
	
inputParser.o: inputParser.cpp ./include/inputParser.h
	$(CPP) $(COMPOPT) $(II) $(LL) -c inputParser.cpp $(LIBS) 

clean:
	rm -rf ./*.o ./libdite2.a
