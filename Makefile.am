SUBDIRS = dynarr dynarr/src


bin_PROGRAMS = primes
primes_SOURCES = test.c primes.c dynarr.h vector.h
primes_CFLAGS = -Idynarr/src/ -Idynarr/vector/src/
primes_LDFLAGS = -static -lm
primes_LDADD = dynarr/src/libdynarr_static.la

