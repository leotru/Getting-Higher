#
# Makefile: Makefile for the NK-landscape implementation.
#
# Wim Hordijk   Last modified: 23 August 2010
#

CC     = g++
CFLAGS = -Wno-deprecated -O3 -c

nk: nk_walk.cpp Random.cpp NK.cpp Random.h NK.h
	$(CC) $(CFLAGS) nk_walk.cpp Random.cpp NK.cpp
	$(CC) -o nk_walk nk_walk.o Random.o NK.o

nkc: nkc_walk.cpp Random.cpp NKC.cpp Random.h NKC.h
	$(CC) $(CFLAGS) nkc_walk.cpp Random.cpp NKC.cpp
	$(CC) -o nkc_walk nkc_walk.o Random.o NKC.o

nkm: nk_measure.cpp Random.cpp NK.cpp Random.h NK.h
	$(CC) $(CFLAGS) nk_measure.cpp Random.cpp NK.cpp
	$(CC) -o nk_measure nk_measure.o Random.o NK.o

nkf: nk_unfold.cpp Random.cpp NK.cpp Random.h NK.h
	$(CC) $(CFLAGS) nk_unfold.cpp Random.cpp NK.cpp
	$(CC) -o nk_unfold nk_unfold.o Random.o NK.o


clean:
	rm -f *~ *.o

#
# EOF: Makefile
#
