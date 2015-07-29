CC = gcc
CXXC = gcc
CFLAGS = -Wall -O3 -g

all: svm-train svm-predict svm-scale

svm-predict: svm-predict.c svm.o
	$(CC) $(CFLAGS) svm-predict.c svm.o -o svm-predict -lm
svm-train: svm-train.c svm.o
	$(CC) $(CFLAGS) svm-train.c svm.o -o svm-train -lm
svm.o: svm.cpp svm.h
	$(CXXC) -c svm.cpp
svm-scale: svm-scale.c
	$(CC) $(CFLAGS) svm-scale.c -o svm-scale
commit:
	LOGNAME=adm cvs commit
clean:
	rm -f *~ svm.o svm-train svm-predict svm-scale
