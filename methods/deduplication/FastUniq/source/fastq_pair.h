/****************************************************************************
 * The 'FASTQ_PAIR' structure group was used to store paired reads and
 * qualities, including basic operation function as well.
 *
 * This file was written by Haibin Xu, December 2011.
 ****************************************************************************/

#include "fastq.h"

#ifndef _FASTQ_PAIR
	typedef struct fastq_pair
		{
			FASTQ_ALL *seq_left; /* the left read */
			FASTQ_ALL *seq_right; /* the right read */
		} FASTQ_PAIR;
	#define _FASTQ_PAIR
#endif

/* create a FASTQ pair. If successful, return the point to it, 
 * otherwise, return NULL. */
FASTQ_PAIR *fastq_pair_create();

/* free the FASTQ pair. If successful, return 0, otherwise return 1. */
int fastq_pair_remove(FASTQ_PAIR *fq_pair);

/* clear the FASTQ pair. If successful, return 0, otherwise return 1. */
int fastq_pair_clear(FASTQ_PAIR *fq_pair);

/* load the left and right reads and qualities for FASTQ pair from file, including description
 * (whether_append_description=1) or not (whether_append_description=0), including quality 
 * (whether_append_quality=1) or not (whether_append_quality=0). If successful, return 0, 
 * otherwise, clear FASTQ pair and return 1. */
int fastq_pair_scanf(FASTQ_PAIR *fq_pair, FILE *fp_left_in, FILE *fp_right_in,
					 int whether_append_description, int whether_append_quality);

/* write the pair-end reads in FASTA or FASTQ format into two output files(format='fa' or 'fq') 
 * or in FASTA format into a single output file(format="fa" and fp_out2==NULL) using the original 
 * description (serial=-1) or the new serial. If successful, return 0, otherwise, return 1. */
int fastq_pair_printf(FASTQ_PAIR *fq_pair, FILE *fp_out1, FILE *fp_out2,
					  char *format, long serial);

/* compare the two FASTQ pairs tightly, if identical, return 0, else if a>b,
 * return 1,  else if a<b, return -1. */
int fastq_pair_compare_tight(FASTQ_PAIR *fq_pair_a, FASTQ_PAIR *fq_pair_b);

/* compare the two FASTQ pairs loosely, if identical, return 0, else if a>b,
 * return 1,  else if a<b, return -1. */
int fastq_pair_compare_loose(FASTQ_PAIR *fq_pair_a, FASTQ_PAIR *fq_pair_b);

/* return the length of left FASTQ sequence in pair, if any error, return -1. */
long fastq_pair_get_left_length(FASTQ_PAIR *fq_pair);

/* return the length of right FASTQ sequence in pair, if any error, return -1. */
long fastq_pair_get_right_length(FASTQ_PAIR *fq_pair);

/* return the length of both left and right FASTQ sequence in pair,
 * if any error, return -1. */
long fastq_pair_get_total_length(FASTQ_PAIR *fq_pair);



