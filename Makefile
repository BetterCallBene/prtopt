CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O0 -g -DNDEBUG
I = -Iinclude -I../sundials/instdir/include

LDFLAGS += -Llib -L../sundials/instdir/lib

LDLIBS += -lm
CS = $(LDFLAGS) -lcsparse -lsundials_cvode -lsundials_nvecserial $(LDLIBS)
CC = icc
MPICC = mpicc
all: clean riccati solver

constraints: src/constraints.c Makefile
	$(CC) $(CF) $(I) -o constraints.o src/constraints.c src/solver.c src/alloc_mem.c src/gen.c src/config.c src/diff.c src/calc_time.c src/nrm.c src/output.c src/rtopt_gen.c src/dyn.c src/cs_ext.c  $(CS)

cost: src/cost.c Makefile
	$(CC) $(CF) $(I) -o cost.o src/cost.c src/solver.c src/alloc_mem.c src/gen.c src/config.c src/diff.c src/calc_time.c src/nrm.c src/output.c src/rtopt_gen.c src/dyn.c src/cs_ext.c  $(CS)

solver: src/solver.c Makefile
	$(CC) $(CF) $(I) -o solver.o src/solver.c src/alloc_mem.c src/gen.c src/config.c src/diff.c src/calc_time.c src/nrm.c src/output.c src/rtopt_gen.c src/dyn.c   $(CS)

dyn: src/dyn.c Makefile
	$(CC) $(CF) $(I) -o dyn.o src/dyn.c src/alloc_mem.c src/gen.c src/config.c src/diff.c src/calc_time.c src/nrm.c src/output.c src/rtopt_gen.c $(CS)

riccati: src/riccati.c Makefile
	$(CC) $(CF) $(I) -o riccati.o src/riccati.c src/solver.c src/alloc_mem.c src/gen.c src/config.c src/diff.c src/calc_time.c src/nrm.c src/output.c src/rtopt_gen.c src/dyn.c src/cs_ext.c  $(CS)

lagrange: src/lagrange.c Makefile
	$(CC) $(CF) $(I) -o lagrange.o src/cost.c src/lagrange.c src/realtimesolver.c src/constraints.c src/riccati.c src/solver.c src/alloc_mem.c src/gen.c src/config.c src/diff.c src/calc_time.c src/nrm.c src/output.c src/rtopt_gen.c src/dyn.c src/cs_ext.c src/test.c  $(CS)	

realtimesolver: src/realtimesolver.c Makefile
	$(CC) $(CF) $(I) -o realtimesolver.o src/cost.c src/lagrange.c src/realtimesolver.c src/constraints.c src/riccati.c src/solver.c src/alloc_mem.c src/gen.c src/config.c src/diff.c src/calc_time.c src/nrm.c src/output.c src/rtopt_gen.c src/dyn.c src/cs_ext.c src/test.c  $(CS)	

realtimesolver_mpi: src/realtimesolver.c Makefile
	$(MPICC) $(CF) $(I) -D WITH_MPI -o realtimesolver.o  src/multipleshooting.c src/cost.c src/lagrange.c src/realtimesolver.c src/constraints.c src/riccati.c src/solver.c src/alloc_mem.c src/gen.c src/config.c src/diff.c src/calc_time.c src/nrm.c src/output.c src/rtopt_gen.c src/dyn.c src/cs_ext.c src/test.c $(CS)	


multi: src/multipleshooting.c Makefile
	$(MPICC) $(CF) $(I) -o multipleshooting.o src/multipleshooting.c src/solver.c src/dyn.c src/rtopt_gen.c src/diff.c src/alloc_mem.c src/gen.c src/config.c src/calc_time.c src/nrm.c src/output.c $(CS)

try: src/try.c Makefile
	$(CC) $(CF) $(I) -D WITH_MPI1 -o try.o src/try.c $(CS)

clean:
	- $(RM) *.o


