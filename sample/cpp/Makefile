CXX = /usr/local/src/gcc-4.9/bin/g++4.9
INCEigen3 = ${HOME}/.local/include
FLAGS = -O0 -std=c++11 -Wl,-rpath=/usr/local/src/gcc-4.9/lib64 -I${INCEigen3} -I../cpp
TARGET = psvd.exe psvd_bin.exe
SRC = psvd.cpp psvd_bin.cpp EmpsvdCore.cpp EmpsvdBinCore.cpp
OBJ = ${SRC:.cpp=.o}

all: ${TARGET}

psvd.exe: psvd.o EmpsvdCore.o
	${CXX} ${FLAGS} -o $@ $^

psvd_bin.exe: psvd_bin.o EmpsvdCore.o EmpsvdBinCore.o
	${CXX} ${FLAGS} -o $@ $^

psvd.o: psvd.cpp EmpsvdCore.o
	${CXX} ${FLAGS} -c $<

psvd_bin.o: psvd_bin.cpp EmpsvdCore.o EmpsvdBinCore.o
	${CXX} ${FLAGS} -c $<

EmpsvdCore.o: EmpsvdCore.cpp
	${CXX} ${FLAGS} -c $<

EmpsvdBinCore.o: EmpsvdBinCore.cpp EmpsvdCore.o
	${CXX} ${FLAGS} -c $<

.PHONY: clean pre

pre:
	cp -v ../data/* ./
	cp -v ../../cpp/* ./

clean:
	${RM} ${TARGET} ${OBJ}

