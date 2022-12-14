include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CFLAGS = -Wall -Werror -g -O0 

.PHONY: clean

all:: clean build

build: 2diff.o chkopts
	echo "******* 2D Diffusion Solver Make begins **********"
	-${CLINKER} ${CFLAGS} -o ../../build/2diff.exe 2diff.o ${PETSC_LIB}
	${RM} 2diff.o
	echo "Make complete: Executable can be found at ../../build/"

clean:: 
	-${RM} edit  ../../build/2diff.exe
	echo "clean completed"


