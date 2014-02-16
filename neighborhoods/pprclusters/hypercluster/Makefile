# Makefile for hypercluster
LIBSPARFUN := ../sparfun
INCLUDE := $(INCLUDE) -I$(LIBSPARFUN)

include ../Makefile.inc

.PHONY : all clean
all: hypercluster

hypercluster.o: hypercluster.h

hypercluster: pagerank.o hypercluster.o $(LIBSPARFUN)/sf_graph.o $(LIBSPARFUN)/sf_load.o \
  $(LIBSPARFUN)/sf_matrix.o $(LIBSPARFUN)/sf_util.o $(LIBSPARFUN)/sf_vec.o
  
clean:
	$(RM) pagerank.o hypercluster.o $(LIBSPARFUN)/sf_graph.o $(LIBSPARFUN)/sf_load.o \
  $(LIBSPARFUN)/sf_matrix.o $(LIBSPARFUN)/sf_util.o $(LIBSPARFUN)/sf_vec.o
