/* This program was used to remove duplicates in paired FASTQ sequences, 
 * which is usually appeared in mate pair libraries.
 *
 * This file and its partner was written by Haibin Xu, December 2011.
 */

#ifndef MAX_FILE_NUMBER
	#define MAX_FILE_NUMBER 1000
#endif

#include <unistd.h>
#include "fastq_pair_array.h"

void fastq_uniq_usage()
{
	fprintf(stderr, "-i : The input file list of paired FSATQ sequence files [FILE IN]\n");
	fprintf(stderr, "        Maximum 1000 pairs\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "     This parameter is used to specify a list of paired sequence files in\n");
	fprintf(stderr, "     FASTQ format as input, in which two adjacent files with reads in the\n");
	fprintf(stderr, "     same order belong to a pair.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-t : Output sequence format [q/f/p]\n");
	fprintf(stderr, "        q : FASTQ format into TWO output files\n");
	fprintf(stderr, "        f : FASTA format into TWO output files\n");
	fprintf(stderr, "        p : FASTA format into ONE output file\n");
	fprintf(stderr, "        default = q\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "     This parameter is used to specify sequence format in output file(s).\n");
	fprintf(stderr, "     FastUniq could output read pairs into two files in either FASTQ [q]\n");
	fprintf(stderr, "     or FASTA [f] format, in which reads in the same order belonging to a\n");
	fprintf(stderr, "     pair. FastUniq could also output read pairs into a single file in\n");
	fprintf(stderr, "     FASTA format [p], in which adjacent reads belonging to a pair.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-o : The first output file [FILE OUT]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-p : The second output file [FILE OUT]\n");
	fprintf(stderr, "     Optional. ONLY required when output sequence format(-t) is specify as\n");
	fprintf(stderr, "     [q] or [f].\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "-c : Types of sequence descriptions for output [0/1]\n");
	fprintf(stderr, "        0 : The raw descriptions\n");
	fprintf(stderr, "        1 : New serial numbers assigned by FastUniq\n");
	fprintf(stderr, "        default = 0\n");
	fprintf(stderr, "\n");
	return;
}

int main (int argc, const char * argv[])
{
	FILE *fp_in_list, *fp_in_left, *fp_in_right, *fp_out_left, *fp_out_right;
	char str_in_left[MAX_FILE_NUMBER][1000], str_in_right[MAX_FILE_NUMBER][1000];
	char str_in_list[1000], str_out_left[1000], str_out_right[1000];
	char s_left[1000], s_right[1000];
	char output_format;
	int description_type;
	int flag_i=0, flag_o=0, flag_t=0, flag_p=0, flag_c=0;
	char ch;
	FASTQ_PAIR *fq_pair;
	FASTQ_PAIR_ARRAY *fq_pair_array, *temp_fq_pair_array;
	long i, seq_pair_count;
	
	if(argc==1)
    {
		fastq_uniq_usage();
		return 1;
	}
	
	/* initializing */
	for(i=0;i<MAX_FILE_NUMBER;i++)
	{
		str_in_left[i][0]='\0';
		str_in_right[i][0]='\0';
	}
	str_in_list[0]='\0';
	str_out_left[0]='\0';
	str_out_right[0]='\0';
	output_format='\0';
	
	/* obtain inputted arguments */
	while((ch=getopt(argc, argv, "i:t:o:p:c:"))!=-1)
    {
        switch(ch)
		{
			case 'i':
                strcpy(str_in_list,optarg);
				if(strcmp(str_in_list,"")!=0)
					flag_i=1;
				else 
				{
					fastq_uniq_usage();
					return 1;
				}
                break;
			case 't':
				if(strlen(optarg)==1)
				{
					if(optarg[0]=='q')
					{
						output_format='q';
						flag_t=1;
						break;
					}
					else if(optarg[0]=='f')
					{
						output_format='f';
						flag_t=1;
						break;
					}
					else if(optarg[0]=='p')
					{
						output_format='p';
						flag_t=1;
						break;
					}					
					else
					{
						fastq_uniq_usage();
						return 1;
					}
				}
				fastq_uniq_usage();
				return 1;
			case 'o':
                strcpy(str_out_left,optarg);
				if(strcmp(str_out_left,"")!=0)
					flag_o=1;
				else 
				{
					fastq_uniq_usage();
					return 1;
				}
                break;
            case 'p':
                strcpy(str_out_right,optarg);
				if(strcmp(str_out_right,"")!=0)
					flag_p=1;
				else 
				{
					fastq_uniq_usage();
					return 1;
				}
                break;
			case 'c':
				if(strlen(optarg)==1)
				{
					if(optarg[0]=='0')
					{
						description_type=0;
						flag_c=1;
						break;
					}
					else if(optarg[0]=='1')
					{
						description_type=1;
						flag_c=1;
						break;
					}
					else
					{
						fastq_uniq_usage();
						return 1;
					}
				}
				fastq_uniq_usage();
				return 1;
            default:
                fastq_uniq_usage();
                break;
		}
	}
	
	/* check inputted arguments */
	if(flag_i==0)
    {
        fprintf(stderr, "Error in input the name of FASTQ file list!\n");
        return 1;
    }
	if(flag_t==0)
		output_format='q';
	if(flag_o==0 || (output_format!='p' && flag_p==0))
	{
		fprintf(stderr, "Error in output sequence file name!\n");
		return 1;
	}
	if(flag_c==0)
		description_type=0;
	
	/* get pair-end FASTQ file list */
	if((fp_in_list=fopen(str_in_list, "r"))==NULL)
    {
        fprintf(stderr, "Error in open FASTQ file list %s for read!\n",
                str_in_list);
        return 1;
    }
	for(i=0; !feof(fp_in_list) && i<MAX_FILE_NUMBER;)
	{
		/* get the file store left FASTQ sequences */
		s_left[0]='\0';
		fgets(s_left, 1000, fp_in_list);
		if(s_left[0]=='\0')
			continue;
		else if(strlen(s_left)>=2 && s_left[strlen(s_left)-1]=='\n')
			s_left[strlen(s_left)-1]='\0';
		else
		{
			fprintf(stderr, "Error in read from FASTQ file list!\n");
			return 1;
		}
		
		/* get the file store right FASTQ sequences */
		s_right[0]='\0';
		fgets(s_right, 1000, fp_in_list);
		if(strlen(s_right)>=2)
		{
			if(s_right[strlen(s_right)-1]=='\n')
				s_right[strlen(s_right)-1]='\0';
		}
		else
		{
			fprintf(stderr, "Error in read from FASTQ file list!\n");
			return 1;
		}
		
		/* append the fiel name to list array */
		strcpy(str_in_left[i], s_left);
		strcpy(str_in_right[i++], s_right);
	}
	fclose(fp_in_list);
	
	/* check the status of pair-end FASTQ files */
	for(i=0;i<MAX_FILE_NUMBER;i++)
	{
		/* check whether list reached the end */
		if(str_in_left[i][0]=='\0')
			break;

		/* check file status */
		if((fp_in_left=fopen(str_in_left[i], "r"))==NULL)
		{
			fprintf(stderr, "Error in open left fastq file %s for read!\n",
					str_in_left[i]);
			return 1;
		}
		fclose(fp_in_left);
		
		if((fp_in_right=fopen(str_in_right[i], "r"))==NULL)
		{
			fprintf(stderr, "Error in open right fastq file %s for read!\n",
					str_in_right[i]);
			return 1;
		}
		fclose(fp_in_right);
	}

	
	/* read all pair-end FASTQ sequences into memory */
	seq_pair_count=0;
	if((fq_pair_array=fastq_pair_array_create())==NULL)
	{
		fprintf(stderr, "Error in allocate enough memory!\n");
		return 1;
	}
	if((temp_fq_pair_array=fastq_pair_array_create())==NULL)
	{
		fprintf(stderr, "Error in allocate enough memory!\n");
		return 1;
	}
	for(i=0;i<MAX_FILE_NUMBER;i++)
	{
		/* check whether list reached the end */
		if(str_in_left[i][0]=='\0')
			break;
		
		/* open inputted pair-end FASTQ file */
		if((fp_in_left=fopen(str_in_left[i], "r"))==NULL)
		{
			fprintf(stderr, "Error in open left fastq file %s for read!\n",
					str_in_left[i]);
			return 1;
		}
		if((fp_in_right=fopen(str_in_right[i], "r"))==NULL)
		{
			fprintf(stderr, "Error in open right fastq file %s for read!\n",
					str_in_right[i]);
			return 1;
		}
		
		/* read sequences */
		for(;!feof(fp_in_left) && !feof(fp_in_right);)
		{
			fq_pair=NULL;
			if((fq_pair=fastq_pair_create())==NULL)
			{
				fprintf(stderr, "Error in allocate enough memory!\n");
				return 1;
			}
			
			if(output_format=='f' || output_format=='p')
			{
				/* NOT require quality */
				if(fastq_pair_scanf(fq_pair, fp_in_left, fp_in_right, description_type==0?1:0, 0)!=0)
				{
					fastq_pair_remove(fq_pair);
					break;
				}
			}
			else
			{
				/* require quality */
				if(fastq_pair_scanf(fq_pair, fp_in_left, fp_in_right, description_type==0?1:0, 1)!=0)
				{
					fastq_pair_remove(fq_pair);
					break;
				}
			}
			
			fastq_pair_array_append(fq_pair, fq_pair_array);
			fastq_pair_array_append(fq_pair, temp_fq_pair_array);
			seq_pair_count++;
		}
		
		if(!feof(fp_in_left) && !feof(fp_in_right))
		{
			fprintf(stderr, "Error in Reading pair-end FASTQ sequence!\n");
			return 1;
		}
	}
	
	/* create memory address index for each BLOCK in a FASTQ_PAIR_ARRAY */
	fastq_pair_array_generate_index(fq_pair_array);
	fastq_pair_array_generate_index(temp_fq_pair_array);
	
	/* sort the pair-end FASTQ sequences */
	fastq_pair_array_sort(fq_pair_array, temp_fq_pair_array, 1, seq_pair_count);
	
	/* open output fastq file */
    if((fp_out_left=fopen(str_out_left, "w"))==NULL)
    {
        fprintf(stderr, "Error in open left fastq file %s for write!\n",
                str_out_left);
        return 1;
    }
	
	if(str_out_right[0]!='\0')
	{
		if((fp_out_right=fopen(str_out_right, "w"))==NULL)
		{
			fprintf(stderr, "Error in open right fastq file %s for write!\n",
					str_out_right);
			return 1;
		}
	}
	
	/* output the sequence in specific format */
	if(output_format=='f')
		fastq_pair_array_printf(fq_pair_array, fp_out_left, fp_out_right, "fa", description_type, 1);
	else if(output_format=='p')
		fastq_pair_array_printf(fq_pair_array, fp_out_left, NULL, "fa", description_type, 1);
	else
		fastq_pair_array_printf(fq_pair_array, fp_out_left, fp_out_right, "fq", description_type, 1);
	
	/* close output files */
	fclose(fp_out_left);
	if(str_out_right[0]!='\0')
		fclose(fp_out_right);
	
//	/* free memory */
//	fastq_pair_array_remove(fq_pair_array);
	
    return 0;
}

