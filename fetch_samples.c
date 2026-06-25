#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <values.h>
#include <time.h>

#define MAXSEQ 100000

typedef union{
	char str[32];
   unsigned long num[8];
}QNAME;

typedef struct{
	QNAME *name;
	QNAME *namedb;
	double sum_q;
	double sum_d;	
	long cluster;
}HIT;


typedef struct{
	char title[256];
	char *seq;
	char *begin;
	char *current;
	int complement;
	int length;
	int coded;
	int index;
   }SEQUENCE;



/*Read a FASTA-format sequence from file*/
int get_next_seq(SEQUENCE *seq, FILE *f){
int i;
unsigned long pos;
char s[2000];
char *buf;
	do{
		if(!fgets(s,1000,f)){
         return(0);
      }
   }while(s[0]!='>');
	strcpy(seq->title, s);
   while(s[strlen(s)-1]!='\n') fgets(s,1000,f);
	buf=(char *)calloc(MAXSEQ, sizeof(char));
	while(1){
		pos=ftell(f);
		if( (!fgets(s,1000,f)) || (s[0]=='>') ||(s[0]=='<')) break;
      s[strlen(s)-1]=0;
      strcat(buf,s);
	}
	fseek(f,pos,SEEK_SET);
   seq->length=strlen(buf);
   seq->seq=seq->begin=seq->current=(char *)calloc(seq->length+10,sizeof(char));
   for(i=0;i<seq->length;i++) seq->seq[i]=toupper(buf[i]);
   free(buf);
   return(1);
}

SEQUENCE *read_seq(FILE *f){
SEQUENCE *seq;
    seq=(SEQUENCE *)calloc(1,sizeof(SEQUENCE));
    get_next_seq(seq,f);
    return(seq);
}

/*Write FASTA-formatted seqience to file*/
void write_seq( SEQUENCE *seq, FILE *f){
int i,j;
	fprintf(f,"%s",seq->title);
   i=j=0;
	while(seq->seq[i++]){
   	fputc(seq->seq[i],f);
      if(j++==60){
      	j=0;
         fprintf(f,"\n");
      }
   }
   fprintf(f,"\n");
}

/*deallocate sequence*/
void free_seq( SEQUENCE *seq){
	if(seq->seq) free(seq->seq);
	free(seq);
}

int startsame(char *str1, char *str2, int range){
int i,flag;
    flag=1;
    for(i=0;i<range;i++){
		if(str1[i]!=str2[i]) flag=0;
    }
    return(flag);
}


void main(int argc, char *argv[]){
unsigned int i,j,num_seq,num_examples,num_clusters,num_hits;
char s1[512], s2[512], s[1024];
HIT **hit_table;
SEQUENCE *seq;
char **example;
time_t t0,t1;
FILE *f,*out;

    if(argc!=4){
	printf("use:\nfetch_samples clusters sequences output\nwhere clusters is a file produced by clustering\n");
	printf("and the sequences are in fasta format (presumably the same file used to BLAST) and output is a file to store results\n");
	exit(-1);
    }

    t0=time(&t0);

    f=fopen(argv[1],"r");
    num_hits=0;
    fgets(s,1000,f);//skip the title line
    while(1){
	if(!fgets(s,1000,f)) break;
	num_hits++;
    }
    rewind(f);
    hit_table=(HIT **)calloc(num_hits,sizeof(HIT *));
    fgets(s,1000,f);//skip the title line
    for(i=0;i<num_hits;i++){
	fgets(s,1000,f);
	hit_table[i]=(HIT *)calloc(1,sizeof(HIT));
	sscanf(s,"%d\t%s\t%s\n",&j,s1,s2);
	s1[31]=0;s2[31]=0;
	hit_table[i]->name=(QNAME *)calloc(1,sizeof(QNAME));
	hit_table[i]->namedb=(QNAME *)calloc(1,sizeof(QNAME));
	strcpy((char *)hit_table[i]->name,s1);
	strcpy((char *)hit_table[i]->namedb,s2);
	hit_table[i]->cluster=j;
    }
    num_clusters=hit_table[i-1]->cluster+1;
    printf("%d hits in %d clusters\n",num_hits, num_clusters);
    fclose(f);

    example=(char **)calloc(num_clusters,sizeof(char *));
    example[0]=(char *)hit_table[0]->name;
    j=1;
    for(i=1;i<num_hits;i++){
	if(hit_table[i]->cluster!=hit_table[i-1]->cluster){
	    example[j++]=(char *)hit_table[i]->name;
	}
    }
    num_examples=j;
    printf("Extracting %d non-redundant examples", num_examples);
    printf("from %s to %s \n", argv[2], argv[3]);


    f=fopen(argv[2],"r");
    out=fopen(argv[3],"w");
    num_seq=0;

    while(1){
	if(!fgets(s,1000,f))break;
	if(s[0]=='>') num_seq++;
    }
    rewind(f);
    printf("from %d sequences\n", num_seq);
    j=0;
    while(1){
	seq=read_seq(f);
	if(!seq->length) break;
	printf("\r%d%% done    ", (int)(j*100)/num_seq);
	j++;
	for(i=0;i<num_examples;i++){
	    if(startsame(seq->title+1,example[i],31)){
		write_seq(seq, out);
		break;
	    }
	}
	free_seq(seq);
    }
    fclose(f);
    fclose(out);
    t1=time(&t1);
    printf("\nFinished, elapsed %ld sec\n",t1-t0);
}
