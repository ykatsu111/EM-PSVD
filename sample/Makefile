
FC = gfortran
FC_FLAGS = -O3 -ffpe-trap=invalid,zero,overflow,underflow -I../fortran
FC_TARGET = sample_f90.exe
FC_SRC = ../fortran/empsvd_static.f90 ../fortran/empsvd_core.f90 sample.f90
FC_OBJ = sample.o empsvd_core.o empsvd_static.o
FC_MOD = sample.mod empsvd_core.mod empsvd_static.mod

CXX = /usr/local/src/gcc-4.9/bin/g++4.9
INCEigen3 = ${HOME}/.local/include
CXX_FLAGS = -O3 -std=c++11 -Wl,-rpath=/usr/local/src/gcc-4.9/lib64 -I${INCEigen3} -I../cpp
CXX_TARGET = sample_c++.exe
CXX_SRC = sample.cpp ../cpp/EmpsvdCore.cpp


all: ${FC_TARGET}

${FC_TARGET}: ${FC_OBJ}
	${FC} ${FC_FLAGS} -o ${FC_TARGET} $^

sample.o: sample.f90 empsvd_core.o empsvd_static.o
	${FC} ${FC_FLAGS} -c $<

empsvd_core.o: ../fortran/empsvd_core.f90 empsvd_static.o
	${FC} ${FC_FLAGS} -c $<

empsvd_static.o: ../fortran/empsvd_static.f90
	${FC} ${FC_FLAGS} -c $<

.PHONY: clean c++

clean:
	${RM} ${CXX_TARGET} ${FC_TARGET} ${FC_OBJ} ${FC_MOD}

c++:
	${CXX} ${CXX_FLAGS} -o ${CXX_TARGET} ${CXX_SRC}
