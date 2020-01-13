CXX = /usr/local/src/gcc-4.9/bin/g++4.9
INCEigen3 = ${HOME}/.local/include
FLAGS = -O3 -std=c++11 -Wl,-rpath=/usr/local/src/gcc-4.9/lib64 -I${INCEigen3}

TARGET = sample.exe
SRC = sample.cpp EmpsvdCore.cpp

all: ${TARGET}

${TARGET}: ${SRC}
	${CXX} ${FLAGS} -o ${TARGET} $^

.PHONY: clean

clean:
	${RM} ${TARGET} ${OBJ}
