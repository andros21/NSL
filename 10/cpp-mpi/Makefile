CC = mpic++
CFLAGS = -Wall -O3 --std=c++11
BOOSTFLAGS = -lboost_system -pthread -lboost_serialization -lboost_mpi

CLASS_PATH = ../../RC
CLASS_O = $(CLASS_PATH)/random.o
CLASS_H = $(CLASS_PATH)/random.h
CLASS_C = $(CLASS_PATH)/random.cpp

TSP_PATH = ../../TSP
TSP_O = $(TSP_PATH)/tsp-mpi.o
TSP_H = $(TSP_PATH)/tsp-mpi.h
TSP_C = $(TSP_PATH)/tsp-mpi.cpp

main.x : main.o $(CLASS_O) $(TSP_O) func.o
	$(CC) $(CLASS_O) $(TSP_O) func.o main.o -o main.x $(BOOSTFLAGS)
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
func.o : func.h func.cpp
	$(CC) -c func.cpp -o func.o $(CFLAGS)
random.o : $(CLASS_H) $(CLASS_C)
	cd $(CLASS_PATH); $(CC) -c random.cpp -o random.o $(CFLAGS)
tsp-mpi.o : $(TSP_H) $(TSP_C)
	cd $(TSP_PATH); $(CC) -c tsp-mpi.cpp -o tsp-mpi.o $(CFLAGS)
clean : clean.local clean.class clean.home clean.tsp-mpi

clean.local :
	rm -f *.o *.gch
clean.class :
	cd $(CLASS_PATH); rm -f *.o *.gch
clean.tsp-mpi :
	cd $(TSP_PATH); rm -f *.o *.gch
clean.home :
	cd ../../; rm -f *.o *.gch
