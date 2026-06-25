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
    unsigned long num[4];
}QNAME;

typedef struct{
    QNAME *name;
    QNAME *namedb;
    int qu;
    int db;
    double sum_qu;
    double sum_db;
    long cluster;
}HIT;

struct DICT{
    QNAME *name;
    double sum;
    struct DICT *bigger;
    struct DICT *smaller;
    struct DICT *same;
    HIT *edge;
};

struct LIST{
    HIT *hit;
    struct LIST *next;
};

struct LIST *insert(struct LIST *list, struct LIST *name){
struct LIST *buf;
    buf=list->next;
    list->next=name;
    name->next=buf;
    return(name);
}

struct LIST *append(struct LIST *list, struct LIST *name){
struct LIST *buf;
    buf=list;
    while(buf->next){
	buf=buf->next;
    }
    buf->next=name;
    return(name);
}

int startsame(char *str1, char *str2, int range){
int i,flag;
    flag=1;
    for(i=0;i<range;i++){
	if(str1[i]!=str2[i]) flag=0;
    }
    return(flag);
}

int same(QNAME *qu, QNAME *db){
int i;
    for(i=0;i<4;i++){
	if(qu->num[i]!=db->num[i]) return(0);
    }
    return(1);
}

static double ruler[4];

void fill_ruler(){
int i,j;
double x;
    x=1.0;
    for(i=0;i<4;i++){
	for(j=0;j<i*64;j++){
	    x*=2.0;
	}
	ruler[i]=x;
    }
}

double sum_up(QNAME *qn){
int i;
double x;
    x=0;
    for(i=1;i<4;i++){
	x+=qn->num[i]*ruler[i];
    }
    return(x);
}

int smaller_of(int x, int y){
    if(x>y) {
	return(y);
    }else{
	return(x);
    }
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
	s[0]=0;
	if(!fgets(s,32,f))break;
	num_taxons++;
    }
    rewind(f);
    taxonlist=(char **)calloc(num_taxons,sizeof(char *));
    for(i=0;i<num_taxons;i++){
	taxonlist[i]=(char *)calloc(32,sizeof(char));
	if(!fgets(taxonlist[i],32,f))break;
	c=strchr(taxonlist[i],'\n');
	if(c) *c=0;
    }
    fclose(f);
    return(num_taxons);
}


//find the ID in the dictionary or the nearest branch to hang it
struct DICT *find( struct DICT *root, struct DICT *name){
    //if already there, return the pointer to the existing structure
    if(root->sum==name->sum){
	return(root);
    }
    //if not, return the pointer to the nearest match (check higher or lower later)
    if(name->sum > root->sum){
	if(!root->bigger){
	    return(root);
	}else{
	    find(root->bigger, name);
	}
    }else{
	if(name->sum < root->sum){
	    if(!root->smaller){
		return(root);
	    }else{
		find(root->smaller, name);
	    }
	}
    }
}

//find the proper branch and hang the new element on it
struct DICT *hang( struct DICT *root, struct DICT *name){
struct DICT *cur;
    cur=find(root,name);
    if(cur->sum==name->sum){
	//propagate cluster marker
	name->edge->cluster=cur->edge->cluster;
	//add one more to the list
	//while(cur->same){
	//	cur=cur->same;
	//}
	name->same=cur->same;
	cur->same=name;
    }else{
	if(name->sum > cur->sum){
	    cur->bigger=name;
	}else{
	    cur->smaller=name;
	}
    }
}


int main(int argc, char *argv[]) {
long i,j,k,flag, num_hits, num_seq, num_clusters,num_sgenes,num_taxons;
int *cluster_size,*code;
FILE *f;
static char *ss, *ss1;
static char s[1000],s1[100],s2[100];
static HIT **hit_table;
static HIT ***cluster;
static struct DICT rootq, rootd;
struct DICT *node, *cur;
struct LIST **clu_hash, *lst;

static time_t t0,t1;

    if(argc!=2){
   		printf("use:\ncl filename \nwhere filename is a tab output of blast\n");
	exit(-1);
    }
    fill_ruler();
    rootq.sum=rootd.sum=ruler[2];

    t0=time(&t0);
    strcpy(s,argv[1]);
    f=fopen(s,"r");
    num_hits=0;
    printf("reading the data\n");


    while(1){
	if(!fgets(s,1000,f)) break;
	num_hits++;
    }
    rewind(f);

    //num_hits=1000000;
    hit_table=(HIT **)calloc(num_hits,sizeof(HIT *));

    for(i=0;i<num_hits;i++){
	if(!fgets(s,1000,f)) break;
	sscanf(s,"%s\t%s\n",s1,s2);
	s1[31]=0;s2[31]=0;
	hit_table[i]=(HIT *)calloc(1,sizeof(HIT));
	printf("%ld%% done     \r",(i*100)/num_hits);
	hit_table[i]->name=(QNAME *)calloc(1,sizeof(QNAME));
	hit_table[i]->namedb=(QNAME *)calloc(1,sizeof(QNAME));
	strcpy((char *)hit_table[i]->name,s1);
	strcpy((char *)hit_table[i]->namedb,s2);
	hit_table[i]->cluster=i;
	//fill up the numeric code
	hit_table[i]->sum_qu=sum_up(hit_table[i]->name);
	hit_table[i]->sum_db=sum_up(hit_table[i]->namedb);
    }
    for(i=1;i<num_hits;i++){
	if(hit_table[i]->sum_qu==hit_table[i-1]->sum_qu) hit_table[i]->cluster=hit_table[i-1]->cluster;
    }

    printf("hanging out left\n");
    for(i=0;i<num_hits;i++){
	node=(struct DICT *)calloc(1,sizeof(struct DICT));
	//hang nodes by query with label propagation
	node->sum=hit_table[i]->sum_qu;
	node->edge=hit_table[i];
	hang(&rootq,node);
	printf("%ld%% done     \r",(i*100)/num_hits);
    }

    printf("hanging out right\n");
    for(i=0;i<num_hits;i++){
	node=(struct DICT *)calloc(1,sizeof(struct DICT));
	//hang nodes by db with label propagation
	node->sum=hit_table[i]->sum_db;
	node->edge=hit_table[i];
	hang(&rootd,node);
	printf("%ld%% done     \r",(i*100)/num_hits);
    }

    printf("connecting dots across\n");
    node=(struct DICT *)calloc(1,sizeof(struct DICT));
    for(i=0;i<num_hits;i++){
	node->sum=hit_table[i]->sum_db;
	node->edge=hit_table[i];
	cur=find(&rootq,node);
	if(cur->sum==node->sum){
	    //propagate the smallest cluster label forward
	    hit_table[i]->cluster=smaller_of(cur->edge->cluster, hit_table[i]->cluster);
        }
	node->sum=hit_table[i]->sum_qu;
	node->edge=hit_table[i];
	cur=find(&rootd,node);
	if(cur->sum==node->sum){
	    hit_table[i]->cluster=smaller_of(cur->edge->cluster, hit_table[i]->cluster);;
        }
	printf("%ld%% done     \r",(i*100)/num_hits);
    }

    free(node);
    fclose(f);

    //build clusters as dynamic lists
    printf("100%% done\nbuilding clusters\n");
    clu_hash=(struct LIST **)calloc(num_hits,sizeof(struct LIST *));
    num_clusters=0;
    flag=-1;
    for(i=0;i<num_hits;i++){
	lst=(struct LIST *)calloc(1,sizeof(struct LIST));
	lst->hit=hit_table[i];
	if(hit_table[i]->cluster>flag){
	    flag=hit_table[i]->cluster;
	    num_clusters++;
	    clu_hash[hit_table[i]->cluster]=lst;
	}else{
	    if(clu_hash[hit_table[i]->cluster]){
		insert(clu_hash[hit_table[i]->cluster],lst);
	    }else{
		clu_hash[hit_table[i]->cluster]=lst;
	    }
	}
	printf("%ld%% done     \r",(i*100)/num_hits);
    }


    printf("100%% done     \n%ld hits\n",num_hits);

    strcpy(s,argv[1]);
    if(strchr(s,'.')) *strchr(s,'.')=0;
    strcat(s,"_cl.txt");

    f=fopen(s,"w");
    printf("writing down the results\n");
    fprintf(f,"cluster#\tfragment_name\tnearest_homolog\n");

    //write down clusters list by list
    k=0;
    for(i=0;i<num_hits;i++){
	if(clu_hash[i]){
	    lst=clu_hash[i];
	    while(1){
		fprintf(f,"%ld\t%s\t%s\n", k, (char *)lst->hit->name, (char *)lst->hit->namedb);
		if(lst->next) {
		    lst=lst->next;
		}else{
		    break;
		}
	    }
	    k++;
	}
	printf("%ld%% done     \r",(i*100)/num_hits);
    }
    num_clusters=k;
/*
    //write down the unsorted list of blast hits with cluster labels
    for(i=0;i<num_hits;i++){
	//if(hit_table[i]->cluster==-1) continue;
	if(hit_table[i]->sum_qu==hit_table[i]->sum_db) continue;
	fprintf(f,"%ld\t%s\t%s\n",hit_table[i]->cluster,(char *)hit_table[i]->name,(char *)hit_table[i]->namedb);
    }
*/
    //write dowwn clusters sorted by number in ascending order
/*    
    printf("Writing the data\n");
    k=0;
    for(i=0;i<num_hits;i++){
	flag=0;
	printf(" %d%% done\r",(i*100)/num_hits);
	for(j=i;j<num_hits;j++){
	    if(hit_table[j]->cluster==i){
		fprintf(f,"%ld\t%s\t%s\n", k, hit_table[j]->name, hit_table[j]->namedb);
		flag++;
		hit_table[j]->cluster=k;
	    }
	}
	if(flag) k++;
    }
*/    

    fclose(f);
    t1=time(&t1);
    printf("cluster labeling done, found %ld records in %ld clusters elapsed %ld sec\n",num_hits, num_clusters, t1-t0);
}


