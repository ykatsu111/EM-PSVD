CXX = /usr/local/src/gcc-4.9/bin/g++4.9
INCEigen3 = ${HOME}/.local/include
CXX_FLAGS = -O3 -std=c++11 -Wl,-rpath=/usr/local/src/gcc-4.9/lib64 -I${INCEigen3} -I../cpp
CXX_TARGET = sample_c++.exe
CXX_SRC = sample.cpp ../cpp/EmpsvdCore.cpp

FC = gfortran
FC_FLAGS = -O3 -ffpe-trap=invalid,zero,overflow,underflow
FC_TARGET = sample_f90.exe
FC_SRC = sample.f90 ../empsvd_core.f90 ../empsvd_static.f90



all: ${CXX_TARGET} ${FC_TARGET}

${CXX_TARGET}: ${CXX_SRC}
	${CXX} ${CXX_FLAGS} -o ${CXX_TARGET} $^

${FC_TARGET}: ${FC_SRC}
	${FC} ${FC_FLAGS} -o ${FC_TARGET} $^

.PHONY: clean

clean:
	${RM} ${CXX_TARGET} ${FC_TARGET}
