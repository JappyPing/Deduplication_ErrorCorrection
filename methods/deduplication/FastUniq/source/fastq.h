/****************************************************************************
 * The 'FASTQ_ALL' structure group was used to store nucleotide sequence in
 * fastq format, including basic operation function as well.
 *
 * This file was written by Haibin Xu, December 2011.
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FASTQ_DESCRIPTION_MAX_LENGTH
	#define FASTQ_DESCRIPTION_MAX_LENGTH 1000
#endif

#ifndef FASTQ_SEQUENCE_MAX_LENGTH
	#define FASTQ_SEQUENCE_MAX_LENGTH 2000
#endif

#ifndef _FASTQ_ALL
	typedef struct fastq_all
		{
			char *description_1; /* the 1st description line */
			char *sequence; /* the sequence*/
			char *description_2; /* the 2nd description line */
			char *quality; /* the quality */
		} FASTQ_ALL;
	#define _FASTQ_ALL
#endif

/* create a FASTQ_ALL sequence. If successful, return the point to it, 
 * otherwise, return NULL. */
FASTQ_ALL *fastq_create();

/* free the FASTQ sequence. If successful, return 0, otherwise return 1. */
int fastq_remove(FASTQ_ALL *fq);

/* clear the FASTQ sequence. If successful, return 0, otherwise return 1. */
int fastq_clear(FASTQ_ALL *fq);

/* get sequence serial from FASTQ description in format '@serial_number'. 
 * If successful return the serial, otherwise return -1. */
long fastq_get_serial(FASTQ_ALL *fq);

/* read a FASTQ sequence from input file, including description (whether_append_description=1)
 * or not (whether_append_description=0), including quality (whether_append_quality=1) or not
 * (whether_append_quality=0). If successful, return 0, otherwise, clear fq
 * and return 1. */
int fastq_scanf(FASTQ_ALL *fq, FILE *fp_in, 
				int whether_append_description, int whether_append_quality);

/* write sequence into output file in FASTQ format(format='fq') or FASTA format(format='fa')
 * using the original description (serial=-1) or the new serial.
 * If successful, return 0, otherwise return 1. */
int fastq_printf(FASTQ_ALL *fq, FILE *fp_out, char *format, long serial);

/* return the length of FASTQ sequence, if any error, return -1. */
long fastq_get_length(FASTQ_ALL *fq);



