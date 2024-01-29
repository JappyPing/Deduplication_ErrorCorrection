/****************************************************************************
 * The 'FASTQ_PAIR_ARRAY' structure group was used to store a array of 
 * paired FASTQ reads, including basic operation function as well.
 *
 * This file was written by Haibin Xu, December 2011.
 ****************************************************************************/

#include "fastq_pair.h"

#ifndef FASTQ_PAIR_ARRAY_BLOCK_SIZE
	#define FASTQ_PAIR_ARRAY_BLOCK_SIZE 100000
#endif

#ifndef _FASTQ_PAIR_ARRAY_BLOCK
	typedef struct fastq_pair_array_block
		{
			FASTQ_PAIR *block[FASTQ_PAIR_ARRAY_BLOCK_SIZE];
			struct fastq_pair_array_block *previous;
			struct fastq_pair_array_block *next;
			long num;
		} FASTQ_PAIR_ARRAY_BLOCK;
	#define _FASTQ_PAIR_ARRAY_BLOCK
#endif

#ifndef _FASTQ_PAIR_ARRAY
	typedef struct fastq_pair_array
		{
			FASTQ_PAIR_ARRAY_BLOCK *array;
			FASTQ_PAIR_ARRAY_BLOCK *last;
			long block_num;
			long fastq_pair_num;
			FASTQ_PAIR_ARRAY_BLOCK **index;
		} FASTQ_PAIR_ARRAY;
	#define _FASTQ_PAIR_ARRAY
#endif

/* create a FASTQ pair array. If successful, return the point to it, 
 * otherwise, return NULL. */
FASTQ_PAIR_ARRAY *fastq_pair_array_create();

/* free the FASTQ pair array. If successful, return 0, otherwise return 1. */
int fastq_pair_array_remove(FASTQ_PAIR_ARRAY *fq_pair_array);

/* append a new FASTQ pair to the array. if successful, return 0, otherwise
 * return 1. */
int fastq_pair_array_append(FASTQ_PAIR *fq_pair,
                            FASTQ_PAIR_ARRAY *fq_pair_array);

/* generate the index for given FASTQ_PAIR, if successful, return 0, otherwise
 * return 1. */
int fastq_pair_array_generate_index(FASTQ_PAIR_ARRAY *fq_pair_array);

/* get double pointer to individual fastq_pair member at specific position
 * in the array, if successful, return the double pointer, otherwise
 * return NULL */
FASTQ_PAIR **fastq_pair_array_get_pointer(FASTQ_PAIR_ARRAY *fq_pair_array,
                                          long position);

/* merge the two sorted part in array, low-middle and middle-high, into a 
 * single sorted order. If successful, return 0, otherwise return. */
int fastq_pair_array_merge(FASTQ_PAIR_ARRAY *fq_pair_array,
                           FASTQ_PAIR_ARRAY *temp_fq_pair_array, 
						   long low, long middle, long high);

/* sort the FASTQ pair array. If successful, return 0, otherwise
 * return 1. */
int fastq_pair_array_sort(FASTQ_PAIR_ARRAY *fq_pair_array, 
                          FASTQ_PAIR_ARRAY *temp_fq_pair_array,
						  long first, long last);

/* write the pair-end reads in the array in FASTA or FASTQ format into two 
 * output files(format='fa' or 'fq')  or in FASTA format into a single output
 * file(format="fa" and fp_out2==NULL) using the original description 
 * (serial_flag=0) or a new serial number(serial_flag=1). Output all sequences
 * (flag_uniq==0), or unique ones(flag_uniq==1). If successful, return 0,
 * otherwise return 1. */
int fastq_pair_array_printf(FASTQ_PAIR_ARRAY *fq_pair_array, FILE *fp_out1,
                            FILE *fp_out2, char *format, int serial_flag, 
                            int flag_uniq);

















