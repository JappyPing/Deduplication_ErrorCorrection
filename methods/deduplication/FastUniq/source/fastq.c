/****************************************************************************
 * The 'FASTQ_ALL' structure group was used to store nucleotide sequence in
 * fastq format, including basic operation function as well.
 *
 * This file was written by Haibin Xu, December 2011.
 ****************************************************************************/

#include "fastq.h"

FASTQ_ALL *fastq_create()
{
	/* create a FASTQ_ALL sequence. If successful, return the point to it, 
	 * otherwise, return NULL/.
	 */
	FASTQ_ALL *fq;
	
	if((fq=(FASTQ_ALL *)malloc(sizeof(FASTQ_ALL)))==NULL)
		return NULL;

	fq->description_1=NULL;
	fq->sequence=NULL;
	fq->description_2=NULL;
	fq->quality=NULL;
	
	return fq;
}

int fastq_remove(FASTQ_ALL *fq)
{
	/* free the FASTQ sequence. If successful, return 0, otherwise return 1.
	 */
	if(fq==NULL)
		return 1;
	
	if(fq->description_1!=NULL)
		free(fq->description_1);
	if(fq->sequence!=NULL)
		free(fq->sequence);
	if(fq->description_2!=NULL)
		free(fq->description_2);
	if(fq->quality!=NULL)
		free(fq->quality);
	
	free(fq);
	
	return 0;
}

int fastq_clear(FASTQ_ALL *fq)
{
    /* clear the FASTQ sequence. If successful, return 0, otherwise return 1.
     */
    if(fq==NULL)
        return 1;
    
	if(fq->description_1!=NULL)
	{
		free(fq->description_1);
		fq->description_1=NULL;
	}
	if(fq->sequence!=NULL)
	{
		free(fq->sequence);
		fq->sequence=NULL;
	}
	if(fq->description_2!=NULL)
	{
		free(fq->description_2);
		fq->description_2=NULL;
	}
	if(fq->quality!=NULL)
	{
		free(fq->quality);
		fq->quality=NULL;
	}
	
    return 0;
}

long fastq_get_serial(FASTQ_ALL *fq)
{
    /* get sequence serial from FASTQ description in format '@serial_number'. 
     * If successful return the serial, otherwise return -1.
     */
    long serial;
    
    if(fq==NULL || fq->description_1==NULL || fq->description_1[0]=='\0')
        return -1;
    
    if((sscanf(fq->description_1, "@%ld", &serial))!=1)
        return -1;
    
    return serial;
}

int fastq_scanf(FASTQ_ALL *fq, FILE *fp_in, 
				int whether_append_description, int whether_append_quality)
{
	/* read a FASTQ sequence from input file, including description (whether_append_description=1)
	 * or not (whether_append_description=0), including quality (whether_append_quality=1) or not
	 * (whether_append_quality=0). If successful, return 0, otherwise, clear fq
	 * and return 1.
	 */
	char description_1[FASTQ_DESCRIPTION_MAX_LENGTH], sequence[FASTQ_SEQUENCE_MAX_LENGTH];
	char description_2[FASTQ_DESCRIPTION_MAX_LENGTH], quality[FASTQ_SEQUENCE_MAX_LENGTH];
	
	char *p_description_1, *p_sequence, *p_description_2, *p_quality;
	
    if(fp_in==NULL || fq==NULL)
        return 1;
    
	fastq_clear(fq);

	/* read the FASTQ sequence */
    fgets(description_1, FASTQ_DESCRIPTION_MAX_LENGTH, fp_in);
    fgets(sequence, FASTQ_SEQUENCE_MAX_LENGTH, fp_in);
    fgets(description_2, FASTQ_DESCRIPTION_MAX_LENGTH, fp_in);
    fgets(quality, FASTQ_SEQUENCE_MAX_LENGTH, fp_in);
	
	/* check whether integrity of the FASTQ sequence */
    if(description_1[0]=='\0' || sequence[0]=='\0' || description_2[0]=='\0' ||
	   quality[0]=='\0' || description_1[0]!='@'|| description_2[0]!='+' ||
	   description_1[strlen(description_1)-1]!='\n' ||
	   sequence[strlen(sequence)-1]!='\n' ||
	   description_2[strlen(description_2)-1]!='\n')
        return 1;

	/* remove return character at the end */
	if(description_1[strlen(description_1)-1]=='\n')
		description_1[strlen(description_1)-1]='\0';
	if(sequence[strlen(sequence)-1]=='\n')
		sequence[strlen(sequence)-1]='\0';
	if(description_2[strlen(description_2)-1]=='\n')
		description_2[strlen(description_2)-1]='\0';
	if(quality[strlen(quality)-1]=='\n')
		quality[strlen(quality)-1]='\0';
	
	/* append the sequence information to fq */
	if((p_sequence=(char *)malloc(strlen(sequence)+1))==NULL)
		return 1;
	strcpy(p_sequence, sequence);
	fq->sequence=p_sequence;
	
	if(whether_append_quality==1)
	{
		if((p_quality=(char *)malloc(strlen(quality)+1))==NULL)
		{
			fastq_clear(fq);
			return 1;
		}
		strcpy(p_quality, quality);
		fq->quality=p_quality;
	}
	if(whether_append_description==1)
	{
		if((p_description_1=(char *)malloc(strlen(description_1)+1))==NULL)
		{
			fastq_clear(fq);
			return 1;
		}
		strcpy(p_description_1, description_1);
		fq->description_1=p_description_1;
		
		if((p_description_2=(char *)malloc(strlen(description_2)+1))==NULL)
		{
			fastq_clear(fq);
			return 1;
		}
		strcpy(p_description_2, description_2);
		fq->description_2=p_description_2;
	}

    return 0;
}

int fastq_printf(FASTQ_ALL *fq, FILE *fp_out, char *format, long serial)
{
    /* write sequence into output file in FASTQ format(format='fq') or FASTA format(format='fa')
	 * using the original description (serial=-1) or the new serial.
	 * If successful, return 0, otherwise return 1.
     */
    if(fp_out==NULL || fq==NULL)
        return 1;
    
	if(strcmp(format, "fq")==0) /* output in FASTQ format */
	{
		if(serial==-1)
		{
			if(fq->description_1!=NULL)
			{
				fputs(fq->description_1, fp_out);
				fputc('\n', fp_out);
				fputs(fq->sequence, fp_out);
				fputc('\n', fp_out);
				fputs(fq->description_2, fp_out);
				fputc('\n', fp_out);
				fputs(fq->quality, fp_out);
				fputc('\n', fp_out);
			}
			else
			{
				fputc('@', fp_out);
				fputc('\n', fp_out);
				fputs(fq->sequence, fp_out);
				fputc('\n', fp_out);
				fputc('+', fp_out);
				fputc('\n', fp_out);
				fputs(fq->quality, fp_out);
				fputc('\n', fp_out);				
			}
		}
		else
		{
			fprintf(fp_out, "@%ld\n",  serial);
			fputs(fq->sequence, fp_out);
			fputc('\n', fp_out);
			fprintf(fp_out, "+%ld\n",  serial);
			fputs(fq->quality, fp_out);
			fputc('\n', fp_out);
		}
	}
	else if(strcmp(format, "fa")==0) /* output in FASTQ format */
	{
		if(serial==-1)
		{
			if(fq->description_1!=NULL)
			{
				fputc('>', fp_out);
				fputs(&(fq->description_1[1]), fp_out);
				fputc('\n', fp_out);
				fputs(fq->sequence, fp_out);
				fputc('\n', fp_out);
			}
			else
			{
				fputc('>', fp_out);
				fputc('\n', fp_out);
				fputs(fq->sequence, fp_out);
				fputc('\n', fp_out);
			}
		}
		else
		{
			fprintf(fp_out, ">%ld\n",  serial);
			fputs(fq->sequence, fp_out);
			fputc('\n', fp_out);
		}
	}
	else
		return 1;
	
    return 0;
}

long fastq_get_length(FASTQ_ALL *fq)
{
	/* return the length of FASTQ sequence, is any error, return -1
	 */
	
	if(fq==NULL)
		return -1;
	if(fq->sequence==NULL)
		return 0;
	return strlen(fq->sequence);
}








