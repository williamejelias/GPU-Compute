# optimised
# CFLAGS = -g -O3 -no-vec -march=native  -D_GNU_SOURCE
# CC = icc

# vectorised
CFLAGS = -g -O3 -qopt-report=5 -qopt-report-phase=vec -march=native  -D_GNU_SOURCE
CC = icc

# CFLAGS = -g -O3 -fopenacc -march=native -D_GNU_SOURCE
# CC = gcc-8.2

# parallelised
# CFLAGS = -g -O3 -fast -acc -Minline -Minfo -ta=multicore -D_GNU_SOURCE
# CC = pgcc

LDFLAGS = -lm

OBJ = optimised-sparsemm.o basic-sparsemm.o utils.o
HEADER = utils.h

.PHONY: clean help check

all: sparsemm

help:
	@echo "Available targets are"
	@echo "  clean: Remove all build artifacts"
	@echo "  check: Perform a simple test of your optimised routines"
	@echo "  sparsemm: Build the sparse matrix-matrix multiplication binary"

clean:
	-rm -f sparsemm $(OBJ)

check: sparsemm
	./sparsemm CHECK

sparsemm: sparsemm.c $(OBJ)
	$(CC) $(CFLAGS) -o $@ $< $(OBJ) $(LDFLAGS)

%.o: %.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $<
