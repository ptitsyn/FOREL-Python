#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <malloc.h>
#include <time.h>


typedef union{
	char str[32];
   unsigned long num[8];
}QNAME;



int smallerof(int x, int y){
	if(x<y){
		return(x);
	}else{
		return(y);
	}
}

typedef struct{
	QNAME *name;
	QNAME *namedb;
	double sum_q;
	double sum_d;	
	long cluster;
}HIT;


double sum(QNAME *s){
int i;
double x;
	for(i=0;i<8;i++){
   	x+=s->num[i];
   }
   return(x);
}

int startsame(char *str1, char *str2, int range){
int i,flag;
    flag=1;
    for(i=0;i<range;i++){
		if(str1[i]!=str2[i]) flag=0;
    }
    return(flag);
}



char **taxonlist;

//Read the list of taxons, return the number of taxons found
//no guardrails: the number and the acronyms of taxons must correspond to the prefixes in the data 
int read_taxonlist(){
int i,num_taxons;
char s[32];
char *c;
FILE *f;
    f=fopen("taxonlist","r");
    if(!f){
	printf("\nmissing the list of taxons\n");
	exit(-1);
    }
    num_taxons=0;
    while(1){
	s[0]=NULL;
	if(!fgets(s,32,f))break;
	num_taxons++;
    }
    rewind(f);
    taxonlist=(char **)calloc(num_taxons,sizeof(char *));
    for(i=0;i<num_taxons;i++){
	taxonlist[i]=(char *)calloc(32,sizeof(char));
	fgets(taxonlist[i],32,f);
	c=strchr(taxonlist[i],'\n');
	if(c) *c=NULL;
    }
    fclose(f);
    return(num_taxons);
}

int tax_code_table[256][256];

void init_tax_code(int num_taxons){
int i;
char a,b;
    for(i=0;i<num_taxons;i++){
	a=taxonlist[i][0];
	b=taxonlist[i][1];
	tax_code_table[a][b]=i;
    }
}

int taxcode(QNAME *name){
char a,b;
    a=name->str[0];
    b=name->str[1];
    return(tax_code_table[a][b]);
}

int main(int argc, char *argv[]) {
unsigned long i,j,k,l,flag, num_hits, num_seq, num_clusters,num_sgenes,num_taxons;
//unsigned long m,mm,gp,st1,st2,st3,st4;
//float id,ev,bs;
int *cluster_size;
FILE *f;
static char *ss, *ss1;
static char s[1000],s1[100],s2[100];
HIT *h;
HIT **buf;
static HIT **hit_table;
static unsigned int **gene_taxon_table;
char **example;
time_t t0,t1;

    if(argc!=2){
   	printf("use:\nsmap_sg filename\nwhere filename is a tab output of cl_sg\n");
   	printf("\nThe list of taxon acronyms must be in file named \"taxonlist\" and each ID in BLAST search output must start with the taxon 2-letter prefix\n");
		exit(-1);
    }

   t0=time(&t0);

    num_taxons=read_taxonlist(/*taxonlist*/);
    init_tax_code(num_taxons);

    strcpy(s,argv[1]);
    f=fopen(s,"r");
    num_hits=0;
    fgets(s,1000,f);//skip the title line
    while(1){
	if(!fgets(s,1000,f)) break;
	num_hits++;
    }
    rewind(f);
    fgets(s,1000,f);//skip the title line

    hit_table=(HIT **)calloc(num_hits,sizeof(HIT *));
    for(i=0;i<num_hits;i++){
	fgets(s,1000,f);
	hit_table[i]=(HIT *)calloc(1,sizeof(HIT));
	sscanf(s,"%ld\t%s\t%s\n",&j,s1,s2);
	s1[31]=0;s2[31]=0;
	hit_table[i]->name=(char *)calloc(1,sizeof(QNAME));
	hit_table[i]->namedb=(char *)calloc(1,sizeof(QNAME));
	strcpy(hit_table[i]->name,s1);
	strcpy(hit_table[i]->namedb,s2);
	hit_table[i]->sum_q=sum(hit_table[i]->name);
	hit_table[i]->sum_d=sum(hit_table[i]->namedb);
	hit_table[i]->cluster=j;
    }
    num_clusters=hit_table[i-1]->cluster+1;
    printf("%d hits in %d clusters\n",num_hits, num_clusters);
    fclose(f);

    cluster_size=(int *)calloc(num_clusters,sizeof(int));
    //allocate taxon table
    gene_taxon_table=(unsigned int **)calloc(num_clusters,sizeof(unsigned int *));
    for(i=0;i<num_clusters;i++) gene_taxon_table[i]=(unsigned int *)calloc(num_taxons,sizeof(unsigned int));
    example=(char **)calloc(num_clusters,sizeof(char *));

    printf("Mapping gene clusters to taxons\n");
    for(i=0;i<num_hits;i++){
	j=hit_table[i]->cluster;
	if(!example[j]) example[j]=hit_table[i]->name;
	k=taxcode(hit_table[i]->name);
	l=taxcode(hit_table[i]->namedb);
	gene_taxon_table[j][k]++;
	if(i==j) continue;
	gene_taxon_table[j][l]++;
	printf("%ld%% done   \r", (i*100)/num_hits);
    }
    printf("100%% done  \n");

    t1=time(&t1);
    printf("ocurrence of %d gene clusters in %d taxons is mapped, elapsed %ld sec\n", num_clusters, num_taxons, t1-t0);
    printf("taxon table ready to print\n");

    strcpy(s,argv[1]);
    if(strchr(s,'.')) *strchr(s,'.')=0;
    strcat(s,"_taxonmap.txt");
    f=fopen(s,"w");
    //output the map, cluster by cluster, where each cluster is the line and columns are taxons
    //each table cell is 0 is cluster/supergene is not found in this taxon or 1 if the supergene is found in the taxon
    fprintf(f,"cluster\t");
    for(i=0;i<num_taxons;i++) fprintf(f,"%s\t",taxonlist[i]);
    fprintf(f,"example\n");
    for(i=0;i<num_clusters;i++){
	if(!example[i]) continue;
	fprintf(f,"sgene%d\t",i+1);
	for(j=0;j<num_taxons;j++){
	    if(gene_taxon_table[i][j]){
		fprintf(f,"1\t");
	    }else{
		fprintf(f,"0\t");
	    }
	}
	fprintf(f,"%s\n",example[i]);
	fprintf(f,"\n");
    }
    fclose(f);
}


