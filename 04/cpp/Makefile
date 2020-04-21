CC = g++
CFLAGS = -Wall -O3 --std=c++11

CLASS_PATH = ../../RC
CLASS_O = $(CLASS_PATH)/random.o
CLASS_H = $(CLASS_PATH)/random.h
CLASS_C = $(CLASS_PATH)/random.cpp

MD_PATH = ../../MD
MD_O = $(MD_PATH)/moldyn.o
MD_H = $(MD_PATH)/moldyn.h
MD_C = $(MD_PATH)/moldyn.cpp

main.x : main.o $(MD_O) func.o
	$(CC) $(MD_O) func.o main.o -o main.x
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
func.o : func.h func.cpp
	$(CC) -c func.cpp -o func.o $(CFLAGS)
random.o : $(CLASS_H) $(CLASS_C)
	cd $(CLASS_PATH); $(CC) -c random.cpp -o random.o $(CFLAGS)
moldyn.o : $(MD_H) $(MD_C)
	cd $(MD_PATH); $(CC) -c moldyn.cpp -o moldyn.o $(CFLAGS)
clean :
	clean.local clean.class clean.moldyn clean.home clean.xyz
clean.local :
	rm -f *.o *.gch
clean.xyz :
	cd ../out/xyz; rm -f *.xyz
clean.class :
	cd $(CLASS_PATH); rm -f *.o *.gch
clean.moldyn :
	cd $(MD_PATH); rm -f *.o *.gch
clean.home :
	cd ../../; rm -f *.o *.gch
