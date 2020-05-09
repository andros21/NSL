CC = g++
CFLAGS = -Wall -O3 --std=c++11

CLASS_PATH = ../../RC
CLASS_O = $(CLASS_PATH)/random.o
CLASS_H = $(CLASS_PATH)/random.h
CLASS_C = $(CLASS_PATH)/random.cpp

I1D_PATH = ../../I1D
I1D_O = $(I1D_PATH)/i1d.o
I1D_H = $(I1D_PATH)/i1d.h
I1D_C = $(I1D_PATH)/i1d.cpp

main.x : main.o $(CLASS_O) $(I1D_O) func.o
	$(CC) $(CLASS_O) $(I1D_O) func.o main.o -o main.x
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
func.o : func.h func.cpp
	$(CC) -c func.cpp -o func.o $(CFLAGS)
i1d.o : $(I1D_H) $(I1D_C)
	cd $(I1D_PATH); $(CC) -c i1d.cpp -o i1d.o $(CFLAGS)
random.o : $(CLASS_H) $(CLASS_C)
	cd $(CLASS_PATH); $(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	clean.local clean.class clean.i1d clean.home
clean.local :
	rm -f *.o *.gch
clean.class :
	cd $(CLASS_PATH); rm -f *.o *.gch
clean.i1d :
	cd $(I1D_PATH); rm -f *.o *.gch
clean.home :
	cd ../../; rm -f *.o *.gch
