#include <stdlib.h>
#include <stdio.h> 
#include <string.h> 
#include <math.h> 
#include <time.h>
#include <R.h>
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define randomize() srand((unsigned)time(NULL));


void ran3(int *idum,int *RN,int *max);
void ran1(int *RN,int *max);
int sort_function(const void *a, const void *b);
 
void distribution(int *idum,int *SAMPLING_CONT,int *LIST_WHOLE_COUNT,
		   int *TOTAL_SAMPLING_COUNT,int *diffcount,int *CHR_NAME,int *START,int *END,
		   int *windowsize);
void distribution1(int *idum,int *SAMPLING_CONT,int *LIST_WHOLE_COUNT,
		   int *TOTAL_SAMPLING_COUNT,int *diffcount,int *CHR_NAME,int *START,int *END,
		   int *windowsize);

void ran1(int *RN,int *max)
{ *RN = rand()%(*max);}

void ran3(int *idum,int *RN, int *max)
{ static int inext, inextp;
  static int ma[56];
  static int iff=0;
  int mj, mk;
  int i,ii,k;
  
  if(*idum <0 || iff==0) {
    iff=1; 
     mj=labs(MSEED-labs(*idum));
     mj %=MBIG;
     ma[55]<-mj;
     for(i=1;i<=54;i++){
       ii=(21*i)%55;
       ma[ii]=mk;
       mk=mj-mk;
       if(mk < MZ) mk += MBIG;
       mj=ma[ii];
     }
     for(k=1;k<=4;k++)
       for(i=1;i<=55;i++) {
	 ma[i] -= ma[1+(i+30)%55];
         if(ma[i] < MZ) ma[i] += MBIG;
       }
     inext=0;
     inextp=31;
     *idum=1;
  }
  if(++inext==56) inext=1;
  if(++inextp ==56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if(mj < MZ) mj += MBIG;
  ma[inext] = mj;
  *RN = mj*FAC*(*max)+1;
}

int sort_function(const void *a, const void *b)
{
        if( *(int *)a < *(int *)b ) return -1;
        else if( *(int *)a > *(int *)b ) return 1;
        else return 0;
}

 
void distribution(int *idum,int *SAMPLING_CONT,int *LIST_WHOLE_COUNT,
		   int *TOTAL_SAMPLING_COUNT,int *diffcount,int *CHR_NAME,int *START,int *END,
                    int *windowsize)
{
  int ii,j,k,temp,temp1,diff,LIST_SAM_COUNT=0,diff_temp;     
  int *ttemp,*LIST_SAM; 
  
  LIST_SAM = Calloc(*SAMPLING_CONT, int); 
  ttemp = Calloc(1, int);  

  for(ii=0; ii<*TOTAL_SAMPLING_COUNT; ii++)
  {   LIST_SAM_COUNT=0;
      for(j=0;j<*SAMPLING_CONT;j++)
      {   ran3(idum,ttemp,LIST_WHOLE_COUNT);
	  for(k=0;k<LIST_SAM_COUNT;k++)
	      if( LIST_SAM[k] == *ttemp ) break;
	  if( k == LIST_SAM_COUNT)
          {   LIST_SAM[LIST_SAM_COUNT] = *ttemp;
	      LIST_SAM_COUNT++;
	  } else {
              j--;
	  }
      }
      qsort((void *)LIST_SAM, LIST_SAM_COUNT, sizeof(LIST_SAM[0]), sort_function); 
      diff_temp=0;    
      for(j=0;j<(*SAMPLING_CONT-1);j++)
      {  temp = LIST_SAM[j];	
	 temp1 = LIST_SAM[j+1];
	 if( CHR_NAME[temp]==CHR_NAME[temp1] )
         {  diff = START[temp1] - END[temp]; 
	    if( diff <*windowsize) 
		diff_temp++;
         }
      } 
      diffcount[ii]=diff_temp;
  }
  Free(LIST_SAM);
  Free(ttemp);
}
    
 
void distribution1(int *idum,int *SAMPLING_CONT,int *LIST_WHOLE_COUNT,
		   int *TOTAL_SAMPLING_COUNT,int *diffcount,int *CHR_NAME,int *START,int *END,
                    int *windowsize)
{
  int ii,j,k,temp,temp1,diff,LIST_SAM_COUNT=0,diff_temp;     
  int *ttemp,*LIST_SAM; 
  
  LIST_SAM = Calloc(*SAMPLING_CONT, int); 
  ttemp = Calloc(1, int);  
  randomize();
  for(ii=0; ii<*TOTAL_SAMPLING_COUNT; ii++)
  {   LIST_SAM_COUNT=0;
      for(j=0;j<*SAMPLING_CONT;j++)
      {   ran1(ttemp,LIST_WHOLE_COUNT);
	  for(k=0;k<LIST_SAM_COUNT;k++)
	      if( LIST_SAM[k] == *ttemp ) break;
	  if( k == LIST_SAM_COUNT)
          {   LIST_SAM[LIST_SAM_COUNT] = *ttemp;
	      LIST_SAM_COUNT++;
	  } else {
              j--;
	  }
      }
      qsort((void *)LIST_SAM, LIST_SAM_COUNT, sizeof(LIST_SAM[0]), sort_function); 
      diff_temp=0;    
      for(j=0;j<(*SAMPLING_CONT-1);j++)
      {  temp = LIST_SAM[j];	
	 temp1 = LIST_SAM[j+1];
	 if( CHR_NAME[temp]==CHR_NAME[temp1] )
         {  diff = START[temp1] - END[temp]; 
	    if( diff <*windowsize) 
		diff_temp++;
         }
      } 
      diffcount[ii]=diff_temp;
  }
  Free(LIST_SAM);
  Free(ttemp);
}
    
