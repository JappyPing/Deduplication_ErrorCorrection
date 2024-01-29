headers= fastq.h fastq_pair.h fastq_pair_array.h
objects= fastq.o fastq_pair.o fastq_pair_array.o

GCC_OPTION= -m64 -O3

all : fastuniq

fastuniq : $(objects) fastq_uniq.o
	gcc -o fastuniq $(objects) fastq_uniq.o $(GCC_OPTION)

fastq.o : fastq.c $(headers)
	gcc -c fastq.c -o fastq.o $(GCC_OPTION)

fastq_pair.o : fastq_pair.c $(headers)
	gcc -c fastq_pair.c -o fastq_pair.o $(GCC_OPTION)

fastq_pair_array.o : fastq_pair_array.c $(headers)
	gcc -c fastq_pair_array.c -o fastq_pair_array.o $(GCC_OPTION)

fastq_uniq.o : fastq_uniq.c $(headers)
	gcc -c fastq_uniq.c -o fastq_uniq.o $(GCC_OPTION)

.PHONY : clean
clean : 
	rm $(objects) fastq_uniq.o
