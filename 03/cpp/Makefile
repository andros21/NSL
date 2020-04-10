CC = g++
CFLAGS = -Wall -O3 --std=c++11

CLASS_PATH = ../../RC
CLASS_O = $(CLASS_PATH)/random.o
CLASS_H = $(CLASS_PATH)/random.h
CLASS_C = $(CLASS_PATH)/random.cpp

main.x : main.o $(CLASS_O) func.o
	$(CC) $(CLASS_O) func.o main.o -o main.x
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
func.o : func.h func.cpp
	$(CC) -c func.cpp -o func.o $(CFLAGS)
random.o : $(CLASS_H) $(CLASS_C)
	cd $(CLASS_PATH); $(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	clean.local clean.class clean.home
clean.local :
	rm -f *.o *.gch
clean.class :
	cd $(CLASS_PATH); rm -f *.o *.gch
clean.home :
	cd ../../; rm -f *.o *.gch
