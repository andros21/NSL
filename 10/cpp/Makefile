CC = g++
CFLAGS = -Wall -O3 --std=c++11

CLASS_PATH = ../../RC
CLASS_O = $(CLASS_PATH)/random.o
CLASS_H = $(CLASS_PATH)/random.h
CLASS_C = $(CLASS_PATH)/random.cpp

TSP_PATH = ../../TSP
TSP_O = $(TSP_PATH)/tsp.o
TSP_H = $(TSP_PATH)/tsp.h
TSP_C = $(TSP_PATH)/tsp.cpp

main.x : main.o $(CLASS_O) $(TSP_O) func.o
	$(CC) $(CLASS_O) $(TSP_O) func.o main.o -o main.x
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
func.o : func.h func.cpp
	$(CC) -c func.cpp -o func.o $(CFLAGS)
random.o : $(CLASS_H) $(CLASS_C)
	cd $(CLASS_PATH); $(CC) -c random.cpp -o random.o $(CFLAGS)
tsp.o : $(TSP_H) $(TSP_C)
	cd $(TSP_PATH); $(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
clean : clean.local clean.class clean.home clean.tsp

clean.local :
	rm -f *.o *.gch
clean.class :
	cd $(CLASS_PATH); rm -f *.o *.gch
clean.tsp :
	cd $(TSP_PATH); rm -f *.o *.gch
clean.home :
	cd ../../; rm -f *.o *.gch
