#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXCLUSTER 256
#define MAXLINES 76043226UL
#define MIDPOINT 5173916648.435934

typedef union{
	char str[32];
   unsigned long num[8];
}QNAME;

typedef struct{
	QNAME query;
   QNAME db;
   float pid;
   int range;
   int mm;
   int gaps;
   int startq;
   int stopq;
   int startd;
   int stopd;
   float eval;
   float bscore;
   char flag;
}BLASTRES;


typedef struct{
	QNAME query;
   QNAME db;
}EDGE;



//EDGE *all_lines[MAXLINES/2];
static QNAME all_names[MAXLINES];

//static unsigned long code[MAXLINES];

BLASTRES *read_blast_line(FILE *f, BLASTRES *r){
static int i,n;
char s[20000],s1[20000],s2[20000];
//BLASTRES *r;
   //r=(BLASTRES *)calloc(1,sizeof(BLASTRES));
  	if(fgets(s,1000,f)){
		sscanf(s,"%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",s1,s2,&r->pid,&r->range,&r->mm,&r->gaps,&r->startq,&r->stopq,&r->startd,&r->stopd,&r->eval,&r->bscore);
		s1[31]=0;s2[31]=0;
      strcpy(&r[i].query,s1);
      strcpy(&r[i].db,s2);
   }else{
   	//free(r);
      return(NULL);
   }
   return(r);
}


int write_blast_line(BLASTRES *r, FILE *f){
static int i;
	fprintf(f,"%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",r->query,r->db,r->pid,r->range,r->mm,r->gaps,r->startq,r->stopq,r->startd,r->stopd,r->eval,r->bscore);
   return(1);
}

int same(QNAME *s1, QNAME *s2){
int i;
	for(i=0;i<8;i++){
   	if(s1->num[i]!=s2->num[i]) return(0);
   }
   return(1);
}

double sum(QNAME *s){
int i;
double x;
    x=0;
    for(i=0;i<8;i++){
	x+=s->num[i];
    }
    return(x);
}

int select_edges( char *fname_in, float MIN_PID, float MIN_E, int MIN_BSCORE, char *fname_out){
static int i,j,n,flag, len, max_len;
double x;
static unsigned long here,total;
static BLASTRES *bl,*b,*bb;
QNAME *qn;
char s[2000];
FILE *in, *out;

    in=fopen(fname_in,"rb");
    out=fopen(fname_out,"w");
    b=(BLASTRES *)calloc(1,sizeof(BLASTRES));
    bl=(BLASTRES *)calloc(1,sizeof(BLASTRES));
    max_len=0;
    bl=read_blast_line(in,bl);
//    printf("%s\t%s %f\t%f\n", bl->query.str, bl->db.str, MIN_PID, bl->pid);
    n=0;
    if(!same(&bl->query,&bl->db)) {
	if(bl->bscore>=MIN_BSCORE){
	    if(bl->eval<=MIN_E){
		if(bl->pid>=MIN_PID) {
		    fprintf(out, "%s\t%s\n", bl->query.str, bl->db.str);
		    n++;
		}
	    }
	}
    }

    while(1){
	b=read_blast_line(in,b);
	if(!b) break;
	if(same(&b->query,&b->db)) continue;
	if(same(&b->query,&bl->query)&&same(&b->db,&bl->db)) {
	    continue;
	}else{
	    bb=bl;bl=b;b=bb;
	    //bl=b;
	}
	if(bl->bscore<MIN_BSCORE) continue;
	if(bl->eval>MIN_E) continue;
	if(bl->pid<MIN_PID) continue;
	fprintf(out, "%s\t%s\n", bl->query.str, bl->db.str);
	n++;
    }

    fclose(in);
    fclose(out);
    free(b);
    free(bl);
    return(n);
}

int main(int argc, char *argv[] ){
int n;
static float MIN_PID;
static int MIN_BSCORE;
static float MIN_E;

    if(argc!=6){
	printf("use:\n\tselect blastres pid eval bitscore outputfile\n\n");
	exit(-1);
    }
    MIN_PID=atof(argv[2]);
    MIN_E=atof(argv[3]);
    MIN_BSCORE=atoi(argv[4]);
    n=select_edges(argv[1], MIN_PID, MIN_E, MIN_BSCORE, argv[5]);
    printf("%d edges selected\n",n);
};





