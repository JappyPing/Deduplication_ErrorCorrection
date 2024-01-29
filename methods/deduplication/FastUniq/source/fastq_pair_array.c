/****************************************************************************
 * The 'FASTQ_PAIR_ARRAY' structure group was used to store a array of 
 * paired FASTQ reads, including basic operation function as well.
 *
 * This file was written by Haibin Xu, December 2011.
 ****************************************************************************/

#include "fastq_pair_array.h"

FASTQ_PAIR_ARRAY *fastq_pair_array_create()
{
    /* create a FASTQ pair array. If successful, return the point to it, 
     * otherwise, return NULL.
     */
    FASTQ_PAIR_ARRAY *fq_pair_array;
    
    if((fq_pair_array=(FASTQ_PAIR_ARRAY *)malloc(sizeof(FASTQ_PAIR_ARRAY)))==NULL)
        return NULL;
    
    if((fq_pair_array->array=
		(FASTQ_PAIR_ARRAY_BLOCK *)malloc(sizeof(FASTQ_PAIR_ARRAY_BLOCK)))==NULL)
	{
		free(fq_pair_array);
		return NULL;
	}
    
	fq_pair_array->last=fq_pair_array->array;
	fq_pair_array->block_num=1;
	fq_pair_array->fastq_pair_num=0;
	
	fq_pair_array->array->previous=NULL;
	fq_pair_array->array->next=NULL;
	fq_pair_array->array->num=0;
	
	fq_pair_array->index=NULL;
	
    return fq_pair_array;
}

int fastq_pair_array_remove(FASTQ_PAIR_ARRAY *fq_pair_array)
{
    /* free the FASTQ pair array. If successful, return 0, otherwise 
     * return 1.
     */
    long i;
	FASTQ_PAIR_ARRAY_BLOCK *fq_pair_array_block;
	
	if(fq_pair_array==NULL)
		return 1;
	
	fq_pair_array_block=fq_pair_array->last;
	for(;fq_pair_array_block!=NULL;)
	{
		for(i=0;i<fq_pair_array_block->num;i++)
			fastq_pair_remove(fq_pair_array_block->block[i]);
		
		fq_pair_array_block=fq_pair_array_block->previous;
	}
	
	if(fq_pair_array->index!=NULL)
		free(fq_pair_array->index);
	
    return 0;
}

int fastq_pair_array_append(FASTQ_PAIR *fq_pair, FASTQ_PAIR_ARRAY *fq_pair_array)
{
    /* append a new FASTQ pair to the array. if successful, return 0, otherwise
     * return 1.
	 */
	FASTQ_PAIR_ARRAY_BLOCK *block_temp;
	
	if(fq_pair_array==NULL || fq_pair==NULL)
		return 1;
	
	if(fq_pair_array->last->num<FASTQ_PAIR_ARRAY_BLOCK_SIZE)
	{
		/* append to the last array_block */
		fq_pair_array->last->block[fq_pair_array->last->num++]=fq_pair;
		fq_pair_array->fastq_pair_num++;
	}
	else
	{
		/* add a new array_block, amd append to it */
		if((block_temp=
			(FASTQ_PAIR_ARRAY_BLOCK *)malloc(sizeof(FASTQ_PAIR_ARRAY_BLOCK)))==NULL)
			return 0;
		
		fq_pair_array->last->next=block_temp;
		block_temp->previous=fq_pair_array->last;
		fq_pair_array->last=block_temp;
		fq_pair_array->block_num++;
		
		block_temp->num=0;
		block_temp->block[block_temp->num++]=fq_pair;
		fq_pair_array->fastq_pair_num++;
	}
	
    return 0;
}

int fastq_pair_array_generate_index(FASTQ_PAIR_ARRAY *fq_pair_array)
{
	/* generate the index for given FASTQ_PAIR, if successful, return 0, otherwise
	 * return 1.
	 */
	FASTQ_PAIR_ARRAY_BLOCK **temp_index;
	FASTQ_PAIR_ARRAY_BLOCK *fq_array_block;
	long i;
	
	if(fq_pair_array==NULL)
		return 1;
	
	if(fq_pair_array->index!=NULL)
	{
		free(fq_pair_array->index);
		fq_pair_array->index=NULL;
	}
	
	if((temp_index=(FASTQ_PAIR_ARRAY_BLOCK **)malloc(sizeof(FASTQ_PAIR_ARRAY_BLOCK *)*(fq_pair_array->block_num)))==NULL)
        return 1;
	
	fq_array_block=fq_pair_array->array;
	for(i=0;i<fq_pair_array->block_num;i++)
	{
		temp_index[i]=fq_array_block;
        fq_array_block=fq_array_block->next;
	}
	
	fq_pair_array->index=temp_index;
	
	return 0;
	
}

FASTQ_PAIR **fastq_pair_array_get_pointer(FASTQ_PAIR_ARRAY *fq_pair_array, long position)
{
    /* get double pointer to individual fastq_pair member at specific position
     * in the array, if successful, return the double pointer, otherwise
     * return NULL
     */
    FASTQ_PAIR_ARRAY_BLOCK *fq_array_block;
    long block_num, num;
    long i;
    
    if(fq_pair_array==NULL || position<=0 || position>fq_pair_array->fastq_pair_num)
        return NULL;
    
    block_num=position/FASTQ_PAIR_ARRAY_BLOCK_SIZE;
    num=position%FASTQ_PAIR_ARRAY_BLOCK_SIZE;
	
    if(num==0)
        num=FASTQ_PAIR_ARRAY_BLOCK_SIZE;
    else
        block_num++;
    
	if(fq_pair_array->index==NULL)
	{
		fq_array_block=fq_pair_array->array;
		for(i=1;i<block_num;i++)
			fq_array_block=fq_array_block->next;
		
		return &fq_array_block->block[num-1];
	}
    else
		return &fq_pair_array->index[block_num-1]->block[num-1];
	
	return NULL;
}

int fastq_pair_array_merge(FASTQ_PAIR_ARRAY *fq_pair_array,
						   FASTQ_PAIR_ARRAY *temp_fq_pair_array, 
						   long low, long middle, long high)
{
    /* merge the two sorted part in array, low-middle and middle-high, into a 
     * single sorted order. If successful, return 0, otherwise return 1.
	 */
    long i, begin1, end1, begin2, end2;
    FASTQ_PAIR **fq_pair_current1, **fq_pair_current2;
    FASTQ_PAIR **temp_fq_pair_current;
	
    if(fq_pair_array==NULL || temp_fq_pair_array==NULL || 
	   low > middle || middle > high || 
	   fq_pair_array->fastq_pair_num!=temp_fq_pair_array->fastq_pair_num)
		return 1;
	
	begin1=low;
    end1=middle;
    begin2=middle+1;
    end2=high;
	
	/* merge processing */
    for(i = low; begin1 <= end1 && begin2 <= end2;i++)
    {
        fq_pair_current1=fastq_pair_array_get_pointer(fq_pair_array, begin1);
        fq_pair_current2=fastq_pair_array_get_pointer(fq_pair_array, begin2);
        
        temp_fq_pair_current=fastq_pair_array_get_pointer(temp_fq_pair_array, i);

        if(fastq_pair_compare_tight(*fq_pair_current1, *fq_pair_current2)<=0)
        {
            *temp_fq_pair_current=*fq_pair_current1;
            begin1++;
        }
        else
        {
            *temp_fq_pair_current=*fq_pair_current2;
            begin2++;
        }
    }
    
	/* moving the remaining data to temp_fq_pair_array */
    if(begin1<=end1)
    {
        for(;begin1<=end1;)
        {
            temp_fq_pair_current=fastq_pair_array_get_pointer(temp_fq_pair_array, i++);
            fq_pair_current1=fastq_pair_array_get_pointer(fq_pair_array, begin1++);
            *temp_fq_pair_current=*fq_pair_current1;
        }
    }
    if(begin2<=end2)
    {
		for(;begin2<=end2;)
		{
			temp_fq_pair_current=fastq_pair_array_get_pointer(temp_fq_pair_array, i++);
			fq_pair_current2=fastq_pair_array_get_pointer(fq_pair_array, begin2++);
			*temp_fq_pair_current=*fq_pair_current2;
		}
    }
    
	/* moving the merged data to original position 'fq_pair_array' */
    for(i=low;i<=high;i++)
    {
        fq_pair_current1=fastq_pair_array_get_pointer(fq_pair_array, i);
        temp_fq_pair_current=fastq_pair_array_get_pointer(temp_fq_pair_array, i);
        *fq_pair_current1=*temp_fq_pair_current;
    }
	
	return 0;
}

int fastq_pair_array_sort(FASTQ_PAIR_ARRAY *fq_pair_array, FASTQ_PAIR_ARRAY *temp_fq_pair_array,
								long first, long last)
{
    /* sort the FASTQ pair array. If successful, return 0, otherwise
     * return 1
     */
	long mid;
    
    if(first<last)
    {
        mid=(first+last)/2;
        fastq_pair_array_sort(fq_pair_array, temp_fq_pair_array, first, mid);
        fastq_pair_array_sort(fq_pair_array, temp_fq_pair_array, mid+1, last);
        fastq_pair_array_merge(fq_pair_array, temp_fq_pair_array, first, mid, last);
    }
    
    return 0;
}

int fastq_pair_array_printf(FASTQ_PAIR_ARRAY *fq_pair_array, FILE *fp_out1, FILE *fp_out2,
                            char *format, int serial_flag, int flag_uniq)
{
    /* write the pair-end reads in the array in FASTA or FASTQ format into two 
     * output files(format='fa' or 'fq')  or in FASTA format into a single output
     * file(format="fa" and fp_out2==NULL) using the original description 
     * (serial_flag=0) or a new serial number(serial_flag=1). Output all sequences
	 * (flag_uniq==0), or unique ones(flag_uniq==1). If successful, return 0,
	 * otherwise return 1.
     */
    long i, k;
    FASTQ_PAIR **temp_fq_pair, **temp_fq_pair_old;
    
    if(flag_uniq==0)
    {
        for(i=1;i<=fq_pair_array->fastq_pair_num;i++)
        {
            temp_fq_pair=fastq_pair_array_get_pointer(fq_pair_array, i);
            
            if(serial_flag==0)
                fastq_pair_printf(*temp_fq_pair, fp_out1, fp_out2, format, -1);
            else
                fastq_pair_printf(*temp_fq_pair, fp_out1, fp_out2, format, i);
        }
    }
    else
    {
		temp_fq_pair_old=fastq_pair_array_get_pointer(fq_pair_array, 1);

		/* the fastq_pair_array contain only one read-pair, output it */
		if(fq_pair_array->fastq_pair_num==1)
		{
			if(serial_flag==0)
				fastq_pair_printf(*temp_fq_pair_old, fp_out1, fp_out2,
								  format, -1);
			else
				fastq_pair_printf(*temp_fq_pair_old, fp_out1, fp_out2,
								  format, k++);
		}
		
		/* compare and output */
        for(i=2, k=1;i<=fq_pair_array->fastq_pair_num;i++)
        {
			temp_fq_pair=fastq_pair_array_get_pointer(fq_pair_array, i);
            if(fastq_pair_compare_loose(*temp_fq_pair_old, *temp_fq_pair)!=0)
            {
                if(serial_flag==0)
                    fastq_pair_printf(*temp_fq_pair_old, fp_out1, fp_out2,
                                      format, -1);
                else
                    fastq_pair_printf(*temp_fq_pair_old, fp_out1, fp_out2,
                                      format, k++);
				
				temp_fq_pair_old=temp_fq_pair;
				
				if(i==fq_pair_array->fastq_pair_num)
				{
					if(serial_flag==0)
						fastq_pair_printf(*temp_fq_pair, fp_out1, fp_out2,
										  format, -1);
					else
						fastq_pair_printf(*temp_fq_pair, fp_out1, fp_out2,
										  format, k++);
				}
            }
            else
            {
                if(fastq_pair_get_left_length(*temp_fq_pair_old) <= fastq_pair_get_left_length(*temp_fq_pair) &&
                   fastq_pair_get_right_length(*temp_fq_pair_old) <= fastq_pair_get_right_length(*temp_fq_pair))
				{
                    temp_fq_pair_old=temp_fq_pair;

					if(i==fq_pair_array->fastq_pair_num)
					{
						if(serial_flag==0)
							fastq_pair_printf(*temp_fq_pair, fp_out1, fp_out2,
											  format, -1);
						else
							fastq_pair_printf(*temp_fq_pair, fp_out1, fp_out2,
											  format, k++);
					}
				}
                else
                {
                    if(serial_flag==0)
                        fastq_pair_printf(*temp_fq_pair_old, fp_out1, fp_out2,
                                          format, -1);
                    else
                        fastq_pair_printf(*temp_fq_pair_old, fp_out1, fp_out2,
                                          format, k++);
                    
                    temp_fq_pair_old=temp_fq_pair;
                }
            }
        }
    }
    return 0;
}



















