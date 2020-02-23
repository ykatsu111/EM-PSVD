FC = gfortran
FLAGS = -O0 -g #-ffpe-trap=invalid,zero,overflow,underflow
TARGET = psvd.exe psvd_bin.exe
SRC = empsvd_static.f90 empsvd_core.f90 empsvd_bin_core.f90 psvd.f90 psvd_bin.f90
OBJ    = ${SRC:.f90=.o}
MOD    = ${SRC:.f90=.mod}

all: ${TARGET}

psvd.exe: psvd.o empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -o $@ $^

psvd_bin.exe: psvd_bin.o empsvd_bin_core.o empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -o $@ $^

psvd.o: psvd.f90 empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -c $<

psvd_bin.o: psvd_bin.f90 empsvd_bin_core.o empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -c $<

empsvd_bin_core.o: empsvd_bin_core.f90 empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -c $<

empsvd_core.o: empsvd_core.f90 empsvd_static.o
	${FC} ${FLAGS} -c $<

empsvd_static.o: empsvd_static.f90
	${FC} ${FLAGS} -c $<


.PHONY: clean pre

pre:
	cp -v ../data/* ./
	cp -v ../../fortran/* ./

clean:
	${RM} ${TARGET} ${OBJ} ${MOD}

