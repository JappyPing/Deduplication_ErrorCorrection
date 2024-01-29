/****************************************************************************
 * The 'FASTQ_PAIR' structure group was used to store paired reads and
 * qualities, including basic operation function as well.
 *
 * This file was written by Haibin Xu, December 2011.
 ****************************************************************************/

#include "fastq_pair.h"

FASTQ_PAIR *fastq_pair_create()
{
	/* create a FASTQ pair. If successful, return the point to it, 
	 * otherwise, return NULL.
	 */
	FASTQ_PAIR *fq_pair;
	
	if((fq_pair=(FASTQ_PAIR *)malloc(sizeof(FASTQ_PAIR)))==NULL)
		return NULL;
	
	fq_pair->seq_left=NULL;
	fq_pair->seq_right=NULL;
	
	return fq_pair;
}

int fastq_pair_remove(FASTQ_PAIR *fq_pair)
{
	/* free the FASTQ pair. If successful, return 0, otherwise return 1.
	 */
	if(fq_pair==NULL)
		return 1;
	
	fastq_pair_clear(fq_pair);
	free(fq_pair);
	
	return 0;
}

int fastq_pair_clear(FASTQ_PAIR *fq_pair)
{
	/* clear the FASTQ pair. If successful, return 0, otherwise return 1.
	 */
	if(fq_pair==NULL)
		return 1;
	
	if(fq_pair->seq_left!=NULL)
	{
		fastq_remove(fq_pair->seq_left);
		fq_pair->seq_left=NULL;
	}
	if(fq_pair->seq_right!=NULL)
	{
		fastq_remove(fq_pair->seq_right);
		fq_pair->seq_right=NULL;
	}
	return 0;
}

int fastq_pair_scanf(FASTQ_PAIR *fq_pair, FILE *fp_left_in, FILE *fp_right_in,
					 int whether_append_description, int whether_append_quality)
{
	/* load the left and right reads and qualities for FASTQ pair from file, including description
	 * (whether_append_description=1) or not (whether_append_description=0), including quality 
	 * (whether_append_quality=1) or not (whether_append_quality=0). If successful, return 0, 
	 * otherwise, clear FASTQ pair and return 1.
	 */
	FASTQ_ALL *fq_left, *fq_right;
	
	if(fq_pair==NULL || fp_left_in==NULL || fp_right_in==NULL)
		return 1;
	
	/* clear the FASTQ_PAIR */
	fastq_pair_clear(fq_pair);
	
	/* create the FASTQ_ALL structure */
	if((fq_left=fastq_create())==NULL)
		return 1;
	if((fq_right=fastq_create())==NULL)
	{
		fastq_remove(fq_left);
		return 1;
	}
	
	if(fastq_scanf(fq_left, fp_left_in, whether_append_description, whether_append_quality)!=0 ||
	   fastq_scanf(fq_right, fp_right_in, whether_append_description, whether_append_quality)!=0)
	{
		fastq_remove(fq_left);
		fastq_remove(fq_right);
		return 1;
	}
	
	fq_pair->seq_left=fq_left;
	fq_pair->seq_right=fq_right;
	
	return 0;
}

int fastq_pair_printf(FASTQ_PAIR *fq_pair, FILE *fp_out1, FILE *fp_out2,
					  char *format, long serial)
{
	/* write the pair-end reads in FASTA or FASTQ format into two output files(format='fa' or 'fq') 
	 * or in FASTA format into a single output file(format="fa" and fp_out2==NULL) using the original 
	 * description (serial=-1) or the new serial. If successful, return 0, otherwise, return 1.
	 */
	if(fq_pair==NULL || fp_out1==NULL)
		return 1;
	
	if((strcmp(format, "fq")==0 && fp_out2!=NULL) ||
		(strcmp(format, "fa")==0 && fp_out2!=NULL))
	{
		fastq_printf(fq_pair->seq_left, fp_out1, format, serial);
		fastq_printf(fq_pair->seq_right, fp_out2, format, serial);
	}
	else if(strcmp(format, "fa")==0 && fp_out2==NULL)
	{
		fastq_printf(fq_pair->seq_left, fp_out1, format, serial);
		fastq_printf(fq_pair->seq_right, fp_out1, format, serial);
	}
	else
		return 1;
	
	return 0;
}

int fastq_pair_compare_tight(FASTQ_PAIR *fq_pair_a, FASTQ_PAIR *fq_pair_b)
{
	/* compare the two FASTQ pairs tightly, if identical, return 0, else if a>b,
	 * return 1,  else if a<b, return -1.
	 */
	char *a_left, *a_right, *b_left, *b_right;
	int i, flag;
	
	/* check whether the sequence read exist */
	if(fq_pair_a==NULL || fq_pair_b==NULL ||
	   fq_pair_a->seq_left==NULL || fq_pair_a->seq_left->sequence==NULL || 
       fq_pair_a->seq_right==NULL || fq_pair_a->seq_right->sequence==NULL || 
       fq_pair_b->seq_left==NULL || fq_pair_b->seq_left->sequence==NULL || 
       fq_pair_b->seq_right==NULL || fq_pair_b->seq_right->sequence==NULL)
		return 1;
    
	/* obtain points to sequence */
	a_left=fq_pair_a->seq_left->sequence;
	a_right=fq_pair_a->seq_right->sequence;
	b_left=fq_pair_b->seq_left->sequence;
	b_right=fq_pair_b->seq_right->sequence;
	
	flag=0;
	for(i=0;;i++)
	{
		if(a_left[i]=='\0' && b_left[i]=='\0')
			break;
		if(a_left[i]==b_left[i])
			continue;
        if(a_left[i]=='\0')
        {
            flag=-1;
            break;
        }
        if(b_left[i]=='\0')
        {
            flag=1;
            break;
        }
		
		switch((int)(a_left[i]>b_left[i]))
		{
			case 1:
				flag=1;
				break;
			case 0:
                flag=-1;
				break;
			default:
				break;
		}
        break;
	}
	
	if(flag==0)
    {
        for(i=0;;i++)
        {
            if(a_right[i]=='\0' && b_right[i]=='\0')
                break;
            if(a_right[i]==b_right[i])
                continue;
            if(a_right[i]=='\0')
            {
                flag=-1;
                break;
            }
            if(b_right[i]=='\0')
            {
                flag=1;
                break;
            }
            
            switch((int)(a_right[i]>b_right[i]))
            {
                case 1:
                    flag=1;
                    break;
                case 0:
                    flag=-1;
                    break;
                default:
                    break;
            }
            break;
        }
    }
	
	return flag;
}

int fastq_pair_compare_loose(FASTQ_PAIR *fq_pair_a, FASTQ_PAIR *fq_pair_b)
{
	/* compare the two FASTQ pairs loosely, if identical, return 0, else if a>b,
	 * return 1,  else if a<b, return -1.
	 */
	char *a_left, *a_right, *b_left, *b_right;
	int i, flag;
	
	/* check whether the sequence read exist */
	if(fq_pair_a==NULL || fq_pair_b==NULL ||
	   fq_pair_a->seq_left==NULL || fq_pair_a->seq_left->sequence==NULL || 
		fq_pair_a->seq_right==NULL || fq_pair_a->seq_right->sequence==NULL || 
		fq_pair_b->seq_left==NULL || fq_pair_b->seq_left->sequence==NULL || 
		fq_pair_b->seq_right==NULL || fq_pair_b->seq_right->sequence==NULL)
		return 1;
		
	/* obtain points to sequence */
	a_left=fq_pair_a->seq_left->sequence;
	a_right=fq_pair_a->seq_right->sequence;
	b_left=fq_pair_b->seq_left->sequence;
	b_right=fq_pair_b->seq_right->sequence;
	
	flag=0;
	for(i=0;;i++)
	{
		if(a_left[i]=='\0' && b_left[i]=='\0')
			break;
		if(a_left[i]==b_left[i])
			continue;
        if(a_left[i]=='\0' || b_left[i]=='\0')
            break;
		
		switch((int)(a_left[i]>b_left[i]))
		{
			case 1:
				flag=1;
				break;
			case 0:
                flag=-1;
				break;
			default:
				break;
		}
        break;
	}
	
	if(flag==0)
    {
        for(i=0;;i++)
        {
            if(a_right[i]=='\0' && b_right[i]=='\0')
                break;
            if(a_right[i]==b_right[i])
                continue;
            if(a_right[i]=='\0' || b_right[i]=='\0')
                break;
            
            switch((int)(a_right[i]>b_right[i]))
            {
                case 1:
                    flag=1;
                    break;
                case 0:
                    flag=-1;
                    break;
                default:
                    break;
            }
            break;
        }
    }
	
	return flag;
}

long fastq_pair_get_left_length(FASTQ_PAIR *fq_pair)
{
	/* return the length of left FASTQ sequence in pair, if any error, return -1.
	 */
	if(fq_pair==NULL)
		return -1;
	return fastq_get_length(fq_pair->seq_left);
}

long fastq_pair_get_right_length(FASTQ_PAIR *fq_pair)
{
	/* return the length of right FASTQ sequence in pair, if any error, return -1.
	 */
	if(fq_pair==NULL)
		return -1;
	return fastq_get_length(fq_pair->seq_right);
}

long fastq_pair_get_total_length(FASTQ_PAIR *fq_pair)
{
	/* return the length of both left and right FASTQ sequence in pair,
	 * if any error, return -1.
	 */
	long left_length, right_length;
	
	if(fq_pair==NULL)
		return -1;
	left_length=fastq_pair_get_left_length(fq_pair);
	right_length=fastq_pair_get_right_length(fq_pair);
	
	if(left_length==-1 || right_length==-1)
		return -1;
	
	return left_length+right_length;
}
