FC = gfortran
FLAGS = -O3 -Wall # -ffpe-trap=invalid,zero,overflow,underflow
TARGET = psvd.exe psvd_bin.exe
SRC = empsvd_static.f90 empsvd_core.f90 psvd.f90 psvd_bin.f90
OBJ    = ${SRC:.f90=.o}
MOD    = ${SRC:.f90=.mod}

all: ${TARGET}

psvd.exe: psvd.o empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -o $@ $^

psvd_bin.exe: psvd_bin.o empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -o $@ $^

psvd.o: psvd.f90 empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -c $<

psvd_bin.o: psvd_bin.f90 empsvd_core.o empsvd_static.o
	${FC} ${FLAGS} -c $<

empsvd_core.o: ./src/empsvd_core.f90 empsvd_static.o
	${FC} ${FLAGS} -c $<

empsvd_static.o: ./src/empsvd_static.f90
	${FC} ${FLAGS} -c $<


.PHONY: clean pre

pre:
	if [ ! -d src ]; then \
		mkdir src; \
	fi
	cp -v ../../fortran/*.f90 ./src

clean:
	${RM} ${TARGET} ${OBJ} ${MOD}
	if [ -d src ]; then \
		rm -rv ./src; \
	fi


