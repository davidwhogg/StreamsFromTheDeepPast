#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

//gsl
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>

/*
  VARIABLE DECLARATIONS
*/
//constants
double halflogtwopi;
//Options:
int K;
char K_opt[]="K";
bool *fixampP, *fixmeanP, *fixcovarP;
char fixamp_opt[]="fixamp";
char fixmean_opt[]="fixmean";
char fixcovar_opt[]="fixcovar";
double avgloglikedata;
double tol;
char tol_opt[]="tol";
long long int maxiter;
char maxiter_opt[]="maxiter";
bool likeonly;
char likeonly_opt[]="likeonly";
double w;
char w_opt[]="w";
int splitnmerge;
char splitnmerge_opt[]="splitnmerge";

//Model parameters:
int d;
int dV;//dim of VV, i.e. d(d+1)/2
struct gaussian{
  double alpha;
  gsl_vector *mm;
  gsl_matrix *VV;
};
struct gaussian *gaussians, *startgaussians;
struct modelbs{
  gsl_vector *bbij;
  gsl_matrix *BBij;
};

//Data parameters:
int N;//number of datapoints
struct datapoint{
  gsl_vector *ww;
  gsl_matrix *SS;
  gsl_matrix *RR;
};
struct datapoint *data, *startdata;

//Log
FILE *logfile;

//EM internals are shared
struct gaussian * newgaussians, * startnewgaussians;
gsl_matrix * qij;
gsl_permutation * p;
gsl_vector * wminusRm, * TinvwminusRm;
gsl_matrix * Tij,* Tij_inv, * VRT, * VRTTinv, * Rtrans, * I;//I is the d dimensional unit matrix
struct modelbs * bs;

//random number generator
gsl_rng * randgen; /* global random number generator */


/*
  FUNCTION DECLARATIONS
*/
bool parse_option(char line[]);
bool read_IC(char line[]);
bool read_data(char line[]);
bool read_till_sep(char curr_value[],FILE *file,char sep);
bool write_model(char outputfilename[]);


/*
  PARSE_OPTION(char line[]): parses a line from the options section of the initial conditions file
*/
bool parse_option(char line[]){
  //Split the option in the name and value
  int ii=0,valuestart;
  char option[12]="", value[1000]="";
  while (line[ii] != '='){
    option[ii]= line[ii];
    ++ii;
  }
  valuestart=++ii;
  while (line[ii] != '\n'){
    value[ii-valuestart]= line[ii];
    ++ii;
  }
  
  //Set parameters
  if (strcmp(option,K_opt) == 0){
    K=atoi(value);
    //Based on this information, allocated some matrices/vectors
    fixampP = (bool *) malloc(K * sizeof (bool));
    fixmeanP = (bool *) malloc(K * sizeof (bool));
    fixcovarP = (bool *) malloc(K * sizeof (bool));
    gaussians = (struct gaussian *) malloc(K * sizeof (struct gaussian) );
    if (fixampP == NULL || fixmeanP == NULL || fixcovarP == NULL || gaussians == NULL){
      printf("Allocation of arrays failed, not enough free memory?\n");
      return -1;
    }
    startgaussians = gaussians;
    
  }
  else if (strcmp(option,maxiter_opt) == 0) maxiter=atoll(value);
  else if (strcmp(option,tol_opt) == 0) tol=atof(value);
  else if (strcmp(option,w_opt) == 0) w=atof(value);
  else if (strcmp(option,likeonly_opt) == 0) likeonly=atoi(value);
  else if (strcmp(option,splitnmerge_opt) == 0) splitnmerge=atoi(value);
  else if (strcmp(option,fixamp_opt) == 0){
    for (ii=0; ii != K; ++ii){
      if (value[ii] == '0') *fixampP=false;
      else *fixampP=true;
      ++fixampP;
    }
    fixampP -= K;
  }
  else if (strcmp(option,fixmean_opt) == 0){
    for (ii=0; ii != K; ++ii){
      if (value[ii] == '0') *fixmeanP=false;
      else *fixmeanP=true;
      ++fixmeanP;
    }
    fixmeanP -= K;
  }
  else if (strcmp(option,fixcovar_opt) == 0){
    for (ii=0; ii != K; ++ii){
      if (value[ii] == '0') *fixcovarP=false;
      else *fixcovarP=true;
      ++fixcovarP;
    }
    fixcovarP -= K;
  }
  else return false;
  
  return true;
}





/*
  READ_IC(): reads the initial conditions file and processes it
 */
bool read_IC(char ICfilename[]){
  
  FILE *ICfile;
  if ( (ICfile= fopen(ICfilename,"r")) == NULL){
    printf ("Opening the initial conditions file failed...\n");
    return false;
  }
  char line[100];
  printf("Reading the options section of the initial conditions file...\n");
  while (fgets(line,100,ICfile) != NULL){
    if (line[0] != '#'){
      if (line[0] == '\n') break;
      if (parse_option(line) == false){
	printf("One of the lines in the options section of the initial conditions file is corrupted\n");
	printf("Please check the initial conditions file and try again\n");
	return false;
      }
    }
  }
  
  printf("Successfully read the options\n");

  printf("Reading the initial model parameters...\n");
  

  //Read first block, establish the dimension d of the modeled quantities
  int countd=0;
  double *mmtemp = (double *) malloc (1000000 * sizeof (double));
  while (fgets(line,100,ICfile) != NULL){
    if (line[0] != '#'){
      if (line[0] == '\n') break;
      *(mmtemp++) = atof(line);
      ++countd;
    }
  }
  //now determine d
  d = (int) (-3 + sqrt(9 + 8 * (countd-1)))/2 ; 
  dV = (int) (d*(d+1)/2);

  //allocate the alpha, mm and VV matrices
  int kk;
  for (kk=0; kk != K; ++kk){
    gaussians->mm = gsl_vector_alloc (d);
    gaussians->VV = gsl_matrix_alloc (d,d);;
    ++gaussians;
  }
  gaussians= startgaussians;
  
  //first map the mmtemp values on the right alpha, mm, VV
  mmtemp -= countd;
  (*gaussians).alpha = *(mmtemp++);
  int dd;
  for (dd=0; dd != d; ++dd)
    gsl_vector_set(gaussians->mm,dd,*(mmtemp++));
  int dd1,dd2;
  for (dd1=0; dd1 != d; ++dd1)
    gsl_matrix_set(gaussians->VV,dd1,dd1,*(mmtemp++));
  for (dd1=0; dd1 != d-1; ++dd1)
    for (dd2=dd1+1; dd2 != d; ++dd2){
      gsl_matrix_set(gaussians->VV,dd1,dd2,*mmtemp);
      gsl_matrix_set(gaussians->VV,dd2,dd1,*mmtemp);
      mmtemp++;
    }

  ++gaussians;
  
  //reallocate mmtemp
  mmtemp -= countd;
  mmtemp = (double *) realloc (mmtemp,countd * sizeof (double) );
  if (mmtemp == NULL){
    printf("Error reallocating memory\n");
    printf("Returning\n");
    return false;
  }

  //Then read the rest of the Gaussians.
  for (kk=1; kk != K; ++kk){
    while (fgets(line,100,ICfile) != NULL){
      if (line[0] != '#'){
	if (line[0] == '\n') break;
	*(mmtemp++) = atof(line);
      }
    }
    mmtemp -=countd;
    (*gaussians).alpha = *(mmtemp++);
    for (dd=0; dd != d; ++dd)
      gsl_vector_set(gaussians->mm,dd,*(mmtemp++));
    for (dd1=0; dd1 != d; ++dd1)
      gsl_matrix_set(gaussians->VV,dd1,dd1,*(mmtemp++));
    for (dd1=0; dd1 != d-1; ++dd1)
      for (dd2=dd1+1; dd2 != d; ++dd2){
	gsl_matrix_set(gaussians->VV,dd1,dd2,*mmtemp);
	gsl_matrix_set(gaussians->VV,dd2,dd1,*mmtemp);
	mmtemp++;
      }
    ++gaussians;
    mmtemp -= countd;
  }
  gaussians = startgaussians;
    
  free(mmtemp);

  fclose(ICfile);

  printf("Successfully read initial model parameters from the initial conditions file\n");


  //Print options
  printf("\nThe options are set to:\n");
  printf("K\t\t=\t");
  printf("%i",K);
  printf("\n");
  printf("maxiter\t\t=\t");
  printf("%lli",maxiter);
  printf("\n");
  printf("tol\t\t=\t");
  printf("%f",tol);
  printf("\n");
  printf("splitnmerge\t=\t");
  printf("%i",splitnmerge);
  printf("\n");
  printf("likeonly\t=\t");
  printf("%i",likeonly);
  printf("\n");
  printf("w\t\t=\t");
  printf("%f",w);
  printf("\n");
  printf("fixamp\t\t=\t");
  int ii;
  for (ii=0; ii != K; ++ii){
    printf("%i",*fixampP);
    if (ii < K-1) printf("\t");
    fixampP++;
  }
  fixampP -= K;
  printf("\n");
  printf("fixmean\t\t=\t");
  for (ii=0; ii != K; ++ii){
    printf("%i",*fixmeanP);
    if (ii < K-1) printf("\t");
    fixmeanP++;
  }
  fixmeanP -= K;
  printf("\n");
  printf("fixcovar\t=\t");
  for (ii=0; ii != K; ++ii){
    printf("%i",*fixcovarP);
    if (ii < K-1) printf("\t");
    fixcovarP++;
  }
  fixcovarP -= K;
  printf("\n");


  //Print the initial model parameters
  printf("\nInitial model parameters used:\n\n");
  for (kk=0; kk != K; ++kk){
    printf("Gaussian ");
    printf("%i",kk);
    printf("\n");
    printf("amp\t=\t");
    printf("%f",(*gaussians).alpha);
    printf("\n");
    printf("mean\t=\t");
    for (dd=0; dd != d; ++dd){
      printf("%f",gsl_vector_get(gaussians->mm,dd));
      if (dd < d-1) printf("\t");
    }
    printf("\n");
    printf("covar\t=\t");
    for (dd1=0; dd1 != d; ++dd1)
      printf("%f\t",gsl_matrix_get(gaussians->VV,dd1,dd1));
    for (dd1=0; dd1 != d-1; ++dd1)
      for (dd2=dd1+1; dd2 != d; ++dd2){
	printf("%f\t",gsl_matrix_get(gaussians->VV,dd1,dd2));
      }
    ++gaussians;
    printf("\n\n");
  }
  gaussians = startgaussians;

  return true;
  
}

/*
  read_till_sep: reads characters into curr_value until sep is reached or until '\n' is reached. If '\n' is the last character, true is returned, otherwise false is returned.
 */
bool read_till_sep(char curr_value[],FILE *file,char sep){
  int vv=0;
  bool not_found_sep = true;
  char curr_char;
  while (not_found_sep){
    curr_char = (char) getc(file);
    if (curr_char == sep) break;
    if (curr_char == '\n') return true;
    curr_value[vv++] = curr_char;
  }

  return false;
}


/*
  READ_DATA(char inputfilename[])
 */
bool read_data(char inputfilename[]){
  
  FILE *inputfile;
  if ( (inputfile= fopen(inputfilename,"r")) == NULL){
    printf ("Opening the data file failed...\n");
    return false;
  }
  
  //First count the number of datapoints
  N=0;
  while (feof(inputfile) ==0) if (getc(inputfile) == '\n') N++;
  
  printf("%i datapoints found in ",N);
  printf(inputfilename);
  printf("\n");
  
  //Then read the actual data
  fseek(inputfile,0,SEEK_SET);

  //allocate data
  data = (struct datapoint *) malloc(N*sizeof (struct datapoint));
  if (data == NULL){
    printf("Allocation of arrays failed, not enough free memory?\n");
    exit(-1);
  } 
  startdata = data;
    
  //current line holds a maximum of d+dV+d*d elements
  //First read the whole line
  double *curr_line = (double *) calloc(d+dV+d*d, sizeof (double) );
  bool end_line;
  int ii=0,dd=0,di;//di is the dimension of the individual datapoints
  while (feof(inputfile) == 0){
    char curr_value[20]="";
    end_line = read_till_sep(curr_value,inputfile,'|');
    *(curr_line++) = atof(curr_value);
    ++dd;
    if (end_line){
      ++ii;
      curr_line -= dd;
      //Determine the dimension of the datapoint, di, from the equation di+di*(di+1)/2+d*di=dd, or, 
      di = (int) ((-(3 + 2 * d) + sqrt((3 + 2 * d)*(3 + 2 * d)+8 * dd))/2);
      //then write data values to memory
      //first allocate the space for the data
      data->ww = gsl_vector_alloc(di);
      data->SS = gsl_matrix_alloc(di,di);
      data->RR = gsl_matrix_alloc(di,d);
      int dd1,dd2;
      for (dd1=0; dd1 != di; ++dd1)
	gsl_vector_set(data->ww,dd1,*(curr_line++));
      for (dd1=0; dd1 != di; ++dd1)
	gsl_matrix_set(data->SS,dd1,dd1,*(curr_line++));
      for (dd1=0; dd1 != di-1; ++dd1)
	for (dd2=dd1+1; dd2 != di; ++dd2){
	  gsl_matrix_set(data->SS,dd1,dd2,*curr_line);
	  gsl_matrix_set(data->SS,dd2,dd1,*curr_line);
	  curr_line++;
	}
      for (dd1=0; dd1 != di; ++dd1)
	for (dd2=0; dd2 != d; ++dd2)
	  gsl_matrix_set(data->RR,dd1,dd2,*(curr_line++));
      curr_line -= dd;
      dd = 0;
      data++;
      if (ii == N) break;
    }
  }
  
  data = startdata;

  fclose(inputfile);
  free(curr_line);

  return true;
}



bool write_model(char outputfilename[]){

  FILE *outputfile;
  if ( (outputfile= fopen(outputfilename,"w+")) == NULL){
    printf ("Opening the data file failed...\n");
    return false;
  }  
  
  //loop over the gaussians and write the results to a file
  int kk,dd,dd1,dd2;
  gaussians = startgaussians;
  for (kk=0; kk != K; ++kk){
    printf("#K=%i\n",kk+1);
    printf("%g\n",(*gaussians).alpha);
    for (dd=0; dd != d; ++dd)
      printf("%g\n",gsl_vector_get(gaussians->mm,dd));
    for (dd1=0; dd1 != d; ++dd1)
      printf("%g\n",gsl_matrix_get(gaussians->VV,dd1,dd1));
    for (dd1=0; dd1 != d-1; ++dd1)
      for (dd2=dd1+1; dd2 != d; ++dd2){
	printf("%g\n",gsl_matrix_get(gaussians->VV,dd1,dd2));
      }
    printf("\n");
    ++gaussians;
  }
  
  gaussians = startgaussians;

  return true;
}

/*
  bovy_det: determinant of matrix
 */
double bovy_det(gsl_matrix * A){
  gsl_permutation * p = gsl_permutation_alloc (A->size1);
  int signum;
  double det;
  gsl_linalg_LU_decomp(A,p,&signum);
  det= gsl_linalg_LU_det(A,signum);
  gsl_permutation_free(p);

  return det;
}




/*
  bovy_isfin : true if finite x
*/
inline bool bovy_isfin(double x){
  return (x > DBL_MAX || x < -DBL_MAX) ? false : true;
}



/*
  minmax : returns the finite minimum and maximum (ie, /NaN in IDL) of a vector that is the row or column of a matrix. Self-explanatory
*/
void minmax(gsl_matrix * q, int row, bool isrow, double * min, double * max){
  *max = -DBL_MAX;
  *min = DBL_MAX;
  int dd;
  double temp;
  if (isrow)
    for (dd = 0; dd != q->size2; ++dd){
      temp  = gsl_matrix_get(q,row,dd);
      if (temp > *max && bovy_isfin(temp))
	*max = temp;
      if (temp < *min && bovy_isfin(temp))
	*min = temp;
    }
  else
    for (dd = 0; dd != q->size1; ++dd){
      temp  = gsl_matrix_get(q,dd,row);
      if (temp > *max && bovy_isfin(temp))
	*max = temp;
      if (temp < *min && bovy_isfin(temp))
	*min = temp;
    }


  return ;
}



/*
  logsum: sum elements that are given as logs, in the row or column of a matrix; Returns answer as a log. Only cares about the ones that make the cut
 
  TO DO: implement such that the whole dynamic range is used
*/
double logsum(gsl_matrix * q, int row, bool isrow){
  double logxmin = log(DBL_MIN);
  double logxmax = log(DBL_MAX);
  int l = (isrow) ? q->size2 : q->size1;

  //First find the maximum and mininum
  double max, min;
  minmax(q,row,isrow,&min,&max);

  //fprintf(logfile,"min=%g\n",min);
  //fprintf(logfile,"max=%g\n",max);


  min *= -1.;
  min += logxmin;
  max *= -1.;
  max += logxmax - log(l);
  (min >  max) ? (max=max) : (max=min);
  double loglike=0.0;
  int dd;
  if (isrow)
    for (dd = 0; dd != q->size2; ++dd)
      loglike += exp(gsl_matrix_get(q,row,dd)+max);
  else
    for (dd = 0; dd != q->size1; ++dd)
      loglike += exp(gsl_matrix_get(q,dd,row)+max);
  
  //fprintf(logfile,"return = %g\n",log(loglike)-max);
  //fprintf(logfile,"max=%g\n",max);

  return log(loglike)-max;
}

/*
  normalize_row(gsl_matrix,int): normalizes the rows of q, these are given as logs!!
*/
double normalize_row(gsl_matrix * q, int row){
  double loglike;
  loglike = logsum(q,row,true);

  int dd;
  for (dd = 0; dd != q->size2; ++dd)
    gsl_matrix_set(q,row,dd,gsl_matrix_get(q,row,dd)-loglike);
  
  return loglike;
}


/*
  EM_step(): does 1 update step in the EM algorithm
*/
bool proj_EM_step(struct datapoint * data, int N, struct gaussian * gaussians, int K,bool * fixamp, bool * fixmean, bool * fixcovar, double * avgloglikedata, bool likeonly, double w){
  *avgloglikedata = 0.0;

  int signum,di;
  double exponent;
  double currqij;

  //Initialize new parameters
  int kk;
  for (kk=0; kk != K; ++kk){
    newgaussians->alpha = 0.0;
    gsl_vector_set_zero(newgaussians->mm);
    gsl_matrix_set_zero(newgaussians->VV);
    ++newgaussians;
  }
  newgaussians= startnewgaussians;
    

  //now loop over data and gaussians to update the model parameters
  int ii, jj;
  for (ii = 0; ii != N; ++ii){
    for (jj = 0; jj != K; ++jj){
      //prepare...
      di = (data->SS)->size1;
      //printf("Datapoint has dimension %i\n",di);
      p = gsl_permutation_alloc (di);
      wminusRm = gsl_vector_alloc (di);
      gsl_vector_memcpy(wminusRm,data->ww);
      TinvwminusRm = gsl_vector_alloc (di);
      Tij = gsl_matrix_alloc(di,di);
      gsl_matrix_memcpy(Tij,data->SS);
      Tij_inv = gsl_matrix_alloc(di,di);
      VRT = gsl_matrix_alloc(d,di);
      VRTTinv = gsl_matrix_alloc(d,di);
      Rtrans = gsl_matrix_alloc(d,di);
      //Calculate Tij
      gsl_matrix_transpose_memcpy(Rtrans,data->RR);
      gsl_blas_dsymm(CblasLeft,CblasUpper,1.0,gaussians->VV,Rtrans,0.0,VRT);//Only the upper right part of VV is calculated
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,data->RR,VRT,1.0,Tij);//This is Tij
      //Calculate LU decomp of Tij and Tij inverse
      gsl_linalg_LU_decomp(Tij,p,&signum);
      gsl_linalg_LU_invert(Tij,p,Tij_inv);
      //Calculate Tijinv*(w-Rm)
      gsl_blas_dgemv(CblasNoTrans,-1.0,data->RR,gaussians->mm,1.0,wminusRm);
      //printf("wminusRm = %f\t%f\n",gsl_vector_get(wminusRm,0),gsl_vector_get(wminusRm,1));
      gsl_blas_dsymv(CblasUpper,1.0,Tij_inv,wminusRm,0.0,TinvwminusRm);
      //printf("TinvwminusRm = %f\t%f\n",gsl_vector_get(TinvwminusRm,0),gsl_vector_get(TinvwminusRm,1));
      gsl_blas_ddot(wminusRm,TinvwminusRm,&exponent);
      //printf("Exponent = %f\nDet = %f\n",exponent,gsl_linalg_LU_det(Tij,signum));
      gsl_matrix_set(qij,ii,jj,log(gaussians->alpha) - di * halflogtwopi - 0.5 * gsl_linalg_LU_lndet(Tij) -0.5 * exponent);//This is actually the log of qij
      //printf("Here we have = %f\n",gsl_matrix_get(qij,ii,jj));
      //Now calculate bij and Bij
      gsl_vector_memcpy(bs->bbij,gaussians->mm);
      gsl_blas_dgemv(CblasNoTrans,1.0,VRT,TinvwminusRm,1.0,bs->bbij);
      //printf("bij = %f\t%f\n",gsl_vector_get(bs->bbij,0),gsl_vector_get(bs->bbij,1));
      gsl_matrix_memcpy(bs->BBij,gaussians->VV);
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,VRT,Tij_inv,0.0,VRTTinv);
      gsl_blas_dgemm(CblasNoTrans,CblasTrans,-1.0,VRTTinv,VRT,1.0,bs->BBij);
      gsl_blas_dsyr(CblasUpper,1.0,bs->bbij,bs->BBij);//This is bijbijT + Bij, which is the relevant quantity
      //Clean up
      gsl_permutation_free (p);
      gsl_vector_free(wminusRm);
      gsl_vector_free(TinvwminusRm);
      gsl_matrix_free(Tij);
      gsl_matrix_free(Tij_inv);
      gsl_matrix_free(VRT);
      gsl_matrix_free(VRTTinv);
      gsl_matrix_free(Rtrans);
      ++gaussians;
      ++bs;
    }
    bs -= K;
    gaussians = startgaussians;
    //Normalize qij properly
    *avgloglikedata += normalize_row(qij,ii);
    if (likeonly){
      bs -= K;
      newgaussians = startnewgaussians;
      ++data;
      continue;
    }
    //printf("qij = %f\t%f\n",gsl_matrix_get(qij,ii,0),gsl_matrix_get(qij,ii,1));
    //printf("avgloglgge = %f\n",*avgloglikedata);
    //Again loop over the gaussians to update the model(can this be more efficient? in any case this is not so bad since K << N)
    for (jj = 0; jj != K; ++jj){
      currqij = exp(gsl_matrix_get(qij,ii,jj));
      //printf("Current qij = %f\n",currqij);
      gsl_vector_scale(bs->bbij,currqij);
      gsl_vector_add(newgaussians->mm,bs->bbij);
      gsl_matrix_scale(bs->BBij,currqij);
      gsl_matrix_add(newgaussians->VV,bs->BBij);
      //printf("bij = %f\t%f\n",gsl_vector_get(bs->bbij,0),gsl_vector_get(bs->bbij,1));
      //printf("Bij = %f\t%f\t%f\n",gsl_matrix_get(bs->BBij,0,0),gsl_matrix_get(bs->BBij,1,1),gsl_matrix_get(bs->BBij,0,1));
      ++newgaussians;
      ++bs;
    }
    bs -= K;
    newgaussians = startnewgaussians;
    ++data;
  }
  data -= N;

  *avgloglikedata /= N;
  if (likeonly)
    return true;

  //Now update the parameters
  //Thus, loop over gaussians again!
  double qj;
  for (jj = 0; jj != K; ++jj){
    qj = exp(logsum(qij,jj,false));//Here we should be careful?
    (qj < DBL_MIN) ? qj = 0: 0;
    //printf("qj = %f\n",qj);
    if (*fixamp != true) {
      gaussians->alpha = qj/N;
      if (qj ==0) {
	*fixamp=1;
	*fixmean=1;
	*fixcovar=1;
	++fixamp;
	++fixmean;
	++fixcovar;
	++gaussians;
	++newgaussians;
	continue;
      }
    }
    if (*fixmean != true){
      gsl_vector_scale(newgaussians->mm,1.0/qj);
      gsl_vector_memcpy(gaussians->mm,newgaussians->mm);
    }
    if (*fixcovar != true){
      if (*fixmean != true)
	gsl_blas_dsyr(CblasUpper,-qj,gaussians->mm,newgaussians->VV);
      else {
	gsl_blas_dsyr(CblasUpper,qj,gaussians->mm,newgaussians->VV);
	gsl_blas_dsyr2(CblasUpper,-qj,gaussians->mm,newgaussians->mm,newgaussians->VV);
      }
      if (w > 0.){
	gsl_matrix_add(newgaussians->VV,I);
      }
      gsl_matrix_scale(newgaussians->VV,1.0/(qj+1.0));
      gsl_matrix_memcpy(gaussians->VV,newgaussians->VV);
    }
    ++fixamp;
    ++fixmean;
    ++fixcovar;
    ++gaussians;
    ++newgaussians;
  }
  newgaussians = startnewgaussians;
  gaussians= startgaussians;
  fixamp -= K;
  fixmean -= K;
  fixcovar -= K;




  return true;
}




/*
  calc_splitnmerge: calculates the splitnmerge hierarchy

  TO DO: - calculate the merge criterium properly
         - figure out what to use as a split criterium
 */
void calc_splitnmerge(int N,struct gaussian * gaussians, int K, gsl_matrix * qij, int snmhierarchy[3][K*(K-1)*(K-2)/2]){
  gsl_matrix * Jmerge = gsl_matrix_alloc(K,K);
  gsl_matrix_set_all(Jmerge,-1.);
  int kk1, kk2, kk, ii,maxsnm= K*(K-1)*(K-2)/2;
  double temp1,temp2,temp;
  for (kk1 = 0; kk1 != K; ++kk1)
    for (kk2 = kk1+1; kk2 != K; ++kk2){
      temp = 0.;
      for (ii=0; ii != N; ++ii){
	//make them all exps
	temp1 = exp(gsl_matrix_get(qij,ii,kk1));
	temp2 = exp(gsl_matrix_get(qij,ii,kk2));
	temp += temp1*temp2;
      }
      gsl_matrix_set(Jmerge,kk1,kk2,temp);
    }
  
  //Then calculate Jsplit
  gsl_vector * Jsplit = gsl_vector_alloc(K);
  gsl_vector * Jsplit_temp = gsl_vector_alloc(K);
  gsl_vector_set_all(Jsplit,-1.);
  for (kk = 0; kk != K; ++kk)
    gsl_vector_set(Jsplit,kk,kk);


  //and put everything in the hierarchy
  size_t maxj, maxk, maxl;
  for (kk1 = 0; kk1 != maxsnm; kk1 += (K-2)){
    gsl_matrix_max_index(Jmerge,&maxj,&maxk);
    gsl_vector_memcpy(Jsplit_temp,Jsplit);
    gsl_vector_set(Jsplit_temp,maxj,-1.);
    gsl_vector_set(Jsplit_temp,maxk,-1.);
    for (kk2=0; kk2 != K-2; ++kk2){
      maxl = gsl_vector_max_index(Jsplit_temp);
      gsl_vector_set(Jsplit_temp,maxl,-1.);
      snmhierarchy[0][kk1+kk2]= maxj;
      snmhierarchy[1][kk1+kk2]= maxk;
      snmhierarchy[2][kk1+kk2]= maxl;
      printf("j = %i, k = %i, l = %i\n",(int)maxj,(int)maxk,(int)maxl);
    }
    //then set it to zero and find the next
    gsl_matrix_set(Jmerge,maxj,maxk,-1.);
  }
    


  //clean up
  gsl_matrix_free(Jmerge);
  gsl_vector_free(Jsplit);
  gsl_vector_free(Jsplit_temp);

  return ;
}



/* 
   bovy_randvec: make d dimensional random vector with typical size
*/
void bovy_randvec(gsl_vector * eps, int d, double length){
  randgen = gsl_rng_alloc(gsl_rng_mt19937);
  length /= sqrt((double)d);
  int dd;
  for (dd = 0; dd != d; ++dd)
    gsl_vector_set(eps,dd,(2.*gsl_rng_uniform(randgen)-1.)*length);

  return;
}



/*
  splitnmerge: splits and merges gaussians
*/
void splitnmergegauss(struct gaussian * gaussians,int K, gsl_matrix * qij, int N, int j, int k, int l){
  //get the gaussians to be split 'n' merged
  d = (gaussians->VV)->size1;//dim of mm
  //j,k,l gaussians
  struct gaussian gaussianj, gaussiank, gaussianl;
  gaussianj.mm = gsl_vector_alloc(d);
  gaussianj.VV = gsl_matrix_alloc(d,d);
  gaussiank.mm = gsl_vector_alloc(d);
  gaussiank.VV = gsl_matrix_alloc(d,d);
  gaussianl.mm = gsl_vector_alloc(d);
  gaussianl.VV = gsl_matrix_alloc(d,d);
  
  gsl_matrix * unitm = gsl_matrix_alloc(d,d);
  gsl_matrix_set_identity(unitm);
  gsl_vector * eps = gsl_vector_alloc(d);
  double qjj,qjk,detVVjl;
  int kk;
  for (kk = 0; kk != K; ++kk){
    if (kk == j){
      gaussianj.alpha = gaussians->alpha;
      gsl_vector_memcpy(gaussianj.mm,gaussians->mm);
      gsl_matrix_memcpy(gaussianj.VV,gaussians->VV);
      qjj = exp(logsum(qij,j,false));
    }
    if (kk == k){
      gaussiank.alpha = gaussians->alpha;
      gsl_vector_memcpy(gaussiank.mm,gaussians->mm);
      gsl_matrix_memcpy(gaussiank.VV,gaussians->VV);
      qjk = exp(logsum(qij,k,false));
    }
    if (kk == l){
      gaussianl.alpha = gaussians->alpha;
      gsl_vector_memcpy(gaussianl.mm,gaussians->mm);
      gsl_matrix_memcpy(gaussianl.VV,gaussians->VV);
    }
    ++gaussians;
  }
  gaussians -= K;

  //merge j & k
  gaussianj.alpha += gaussiank.alpha;
  if (qjk == 0. && qjj == 0){
    gsl_vector_add(gaussianj.mm,gaussiank.mm);
    gsl_vector_scale(gaussianj.mm,0.5);
    gsl_matrix_add(gaussianj.VV,gaussiank.VV);
    gsl_matrix_scale(gaussianj.VV,0.5);
  }
  else{
    gsl_vector_scale(gaussianj.mm,qjj/(qjj+qjk));
    gsl_vector_scale(gaussiank.mm,qjk/(qjj+qjk));
    gsl_vector_add(gaussianj.mm,gaussiank.mm);
    gsl_matrix_scale(gaussianj.VV,qjj/(qjj+qjk));
    gsl_matrix_scale(gaussiank.VV,qjk/(qjj+qjk));
    gsl_matrix_add(gaussianj.VV,gaussiank.VV);
  }

  //split l
  gaussianl.alpha /= 2.;
  gaussiank.alpha = gaussianl.alpha;
  detVVjl = bovy_det(gaussianl.VV);
  detVVjl= pow(detVVjl,1./d);
  gsl_matrix_scale(unitm,detVVjl);
  gsl_matrix_memcpy(gaussiank.VV,unitm);
  gsl_matrix_memcpy(gaussianl.VV,unitm);
  gsl_vector_memcpy(gaussiank.mm,gaussianl.mm);
  bovy_randvec(eps,d,detVVjl);
  gsl_vector_add(gaussiank.mm,eps);
  bovy_randvec(eps,d,detVVjl);
  gsl_vector_add(gaussianl.mm,eps);
  
  //copy everything back into the right gaussians
  for (kk = 0; kk != K; ++kk){
    if (kk == j){
      gaussians->alpha = gaussianj.alpha;
      gsl_vector_memcpy(gaussians->mm,gaussianj.mm);
      gsl_matrix_memcpy(gaussians->VV,gaussianj.VV);
    }
    if (kk == k){
      gaussians->alpha = gaussiank.alpha;
      gsl_vector_memcpy(gaussians->mm,gaussiank.mm);
      gsl_matrix_memcpy(gaussians->VV,gaussiank.VV);
    }
    if (kk == l){
      gaussians->alpha = gaussianl.alpha;
      gsl_vector_memcpy(gaussians->mm,gaussianl.mm);
      gsl_matrix_memcpy(gaussians->VV,gaussianl.VV);
    }
    ++gaussians;
  }
  gaussians -= K;

  //cleanup
  gsl_matrix_free(unitm);
  gsl_vector_free(eps);


  return ;
}



/*
  EM: Expectation-maximization steps
  cals EM_step & ...
*/
bool proj_EM(struct datapoint * data, int N, struct gaussian * gaussians, int K,bool * fixamp, bool * fixmean, bool * fixcovar, double * avgloglikedata, double tol,long long int maxiter, bool likeonly, double w, int splitnmerge){


  double diff = 2 * tol, oldavgloglikedata;
  int niter = 0;
  halflogtwopi  = 0.5 * log(8. * atan(1.0));
  while ( diff > tol && niter < maxiter){
    proj_EM_step(data,N,gaussians,K,fixamp,fixmean,fixcovar,avgloglikedata,likeonly,w);
      printf("Got here, loglike = %f\n",*avgloglikedata);
    if (niter > 0) diff = *avgloglikedata - oldavgloglikedata;
    oldavgloglikedata = *avgloglikedata;
    if (likeonly) break;
    ++niter;
    //write_model("result.dat");
  }
  

  return true;
}




/*
  EMwsnm: Outer loop, runs EM, computes splitnmerge, splits and merges, then runs partial EM followed by full EM. And so on
*/
void proj_gauss_mixtures(struct datapoint * data, int N, struct gaussian * gaussians, int K,bool * fixamp, bool * fixmean, bool * fixcovar, double * avgloglikedata, double tol,long long int maxiter, bool likeonly, double w, int splitnmerge){
  //Allocate some memory
  startgaussians = gaussians;
  d = (gaussians->VV)->size1;//dim of mm
  //allocate the newalpha, newmm and newVV matrices
  newgaussians = (struct gaussian *) malloc(K * sizeof (struct gaussian) );
  startnewgaussians = newgaussians;
  int kk;
  for (kk=0; kk != K; ++kk){
    newgaussians->alpha = 0.0;
    newgaussians->mm = gsl_vector_calloc (d);
    newgaussians->VV = gsl_matrix_calloc (d,d);
    ++newgaussians;
  }
  newgaussians= startnewgaussians;
  //allocate the q_ij matrix
  qij = gsl_matrix_alloc(N,K);;
  I = gsl_matrix_alloc(d,d);
  gsl_matrix_set_identity(I);//Unit matrix
  gsl_matrix_scale(I,w);//scaled to w
  //Also take care of the bbij's and the BBij's
  bs = (struct modelbs *) malloc(K * sizeof (struct modelbs) );
  for (kk = 0; kk != K; ++kk){
    bs->bbij = gsl_vector_alloc (d);
    bs->BBij = gsl_matrix_alloc (d,d);
    ++bs;
  }
  bs -= K;

  double EMavgloglikedata;

  proj_EM(data,N,gaussians,K,fixamp,fixmean,fixcovar,avgloglikedata,tol,maxiter,likeonly,w,splitnmerge);


  //Run splitnmerge
  struct gaussian * oldgaussians;
  bool weretrying = true;
  if (likeonly || splitnmerge == 0)
    ;
  else {
    int maxsnm = K*(K-1)*(K-1)/2;
    int snmhierarchy[3][maxsnm];
    while (weretrying){
      //store avgloglike from normal EM and model parameters
      EMavgloglikedata = *avgloglikedata;
      oldgaussians = (struct gaussian *) malloc(K * sizeof (struct gaussian) );
      for (kk=0; kk != K; ++kk){
	oldgaussians->alpha = gaussians->alpha;
	oldgaussians->mm = gsl_vector_calloc (d);
	gsl_vector_memcpy(oldgaussians->mm,gaussians->mm);
	oldgaussians->VV = gsl_matrix_calloc (d,d);
	gsl_matrix_memcpy(oldgaussians->VV,gaussians->VV);
	++oldgaussians;
	++gaussians;
      }
      gaussians -= K;
      oldgaussians -= K;
      //Then calculate the splitnmerge hierarchy
      calc_splitnmerge(N,gaussians,K,qij,snmhierarchy);
      //Then go through this hierarchy
      kk=0;
      while (kk != splitnmerge && kk != maxsnm){
	//splitnmerge
	splitnmergegauss(gaussians,K,qij,N,snmhierarchy[0][kk],snmhierarchy[1][kk],snmhierarchy[2][kk]);
	//partial EM
	//Full EM
	//Better?
	++kk;
      }
      weretrying = false;
    }

    
    
  }


  
  //Free memory
  gsl_matrix_free(I);
  gsl_matrix_free(qij);
  for (kk = 0; kk != K; ++kk){
    gsl_vector_free(bs->bbij);
    gsl_matrix_free(bs->BBij);
    ++bs;
  }
  bs -= K;
  free(bs);
  for (kk=0; kk != K; ++kk){
    gsl_vector_free(newgaussians->mm);
    gsl_matrix_free(newgaussians->VV);
    ++newgaussians;
  }
  newgaussians= startnewgaussians;
  free(newgaussians);
  gsl_rng_free(randgen);

}






/*
  EM_IDL : takes input from IDL, allocates everything appropriately and runs EM then updates the gaussians and stuff in IDL and returns
  
  return int: code of success



TO DO:

- write log messages to the logfile
- Implement the split and merge algorithm: write function which computes the list, which gets called from EM
- return avgloglike's after each EM-step? For convergence showing purposes, put it in the log
- That's about it.
- Oh right, TEST EXTENSIVELY, make sure covariances stay symmetrical (might not be garantueed by the gsl_... I'm using)
- Start splitting up this one file into several files, i.e:
   - put the algorithm in a separate file, write a header
   - put the IDL function in it's own file
   - put the input/output in an in/output file
   - main should only contain main and several headers
   - That's about it?

 */
int EM_IDL(double * ydata, double * ycovar, double * projection, int N, int dy, double * amp, double * xmean, double * xcovar, int d, int K, char * fixamp, char * fixmean, char * fixcovar, double * avgloglikedata, double tol, int maxiter, char likeonly, double w, char * logfilename, int slen, int splitnmerge){
  bool keeplog = true;
  char logname[slen+1];
  int ss;
  if (*logfilename == 0 || likeonly != 0 || slen == 0)
    keeplog = false;
  else {
    for (ss = 0; ss != slen; ++ss)
      logname[ss] = (char) *(logfilename++);
    logfilename -= slen;
    logname[slen] = '\0';
  }

  //  return logfilename;
  
  if (keeplog) {
    logfile = fopen(logname,"a");
    if (logfile == NULL) return -1;
  }



  if (keeplog){
    time_t now;
    time(&now);
    fprintf(logfile,"----------------------------------\n");
    fprintf(logfile,"\n%s\n",asctime(localtime(&now)));
    fprintf(logfile,"----------------------------------\n");
    fflush(logfile);
  }
  
  //Copy everything into the right formats
  struct datapoint * data = (struct datapoint *) malloc( N * sizeof (struct datapoint) );
  struct gaussian * gaussians = (struct gaussian *) malloc (K * sizeof (struct gaussian) );

  
  int ii, jj,dd1,dd2;
  for (ii = 0; ii != N; ++ii){
    data->ww = gsl_vector_alloc(dy);
    data->SS = gsl_matrix_alloc(dy,dy);
    data->RR = gsl_matrix_alloc(dy,d);
    for (dd1 = 0; dd1 != dy; ++dd1)
      gsl_vector_set(data->ww,dd1,*(ydata++));
    for (dd1 = 0; dd1 != dy; ++dd1)
      for (dd2 = 0; dd2 != dy; ++dd2)
	gsl_matrix_set(data->SS,dd1,dd2,*(ycovar++));
    for (dd1 = 0; dd1 != dy; ++dd1)
      for (dd2 = 0; dd2 != d; ++dd2)
	gsl_matrix_set(data->RR,dd1,dd2,*(projection++));
    ++data;
  }
  data -= N;
  ydata -= N*dy;
  ycovar -= N*dy*dy;
  projection -= N*dy*d;


  for (jj = 0; jj != K; ++jj){
    gaussians->mm = gsl_vector_alloc(d);
    gaussians->VV = gsl_matrix_alloc(d,d);
    gaussians->alpha = *(amp++);
    for (dd1 = 0; dd1 != d; ++dd1)
      gsl_vector_set(gaussians->mm,dd1,*(xmean++));
    for (dd1 = 0; dd1 != d; ++dd1)
      for (dd2 = 0; dd2 != d; ++dd2)
	gsl_matrix_set(gaussians->VV,dd1,dd2,*(xcovar++));
    ++gaussians;
  }
  gaussians -= K;
  amp -= K;
  xmean -= K*d;
  xcovar -= K*d*d;




  //Print the initial model parameters to the logfile
  int kk;
  if (keeplog){
    fprintf(logfile,"\nInitial model parameters used:\n\n");
    for (kk=0; kk != K; ++kk){
      fprintf(logfile,"Gaussian ");
      fprintf(logfile,"%i",kk);
      fprintf(logfile,"\n");
      fprintf(logfile,"amp\t=\t");
      fprintf(logfile,"%f",(*gaussians).alpha);
      fprintf(logfile,"\n");
      fprintf(logfile,"mean\t=\t");
      for (dd1=0; dd1 != d; ++dd1){
	fprintf(logfile,"%f",gsl_vector_get(gaussians->mm,dd1));
	if (dd1 < d-1) fprintf(logfile,"\t");
      }
      fprintf(logfile,"\n");
      fprintf(logfile,"covar\t=\t");
      for (dd1=0; dd1 != d; ++dd1)
	fprintf(logfile,"%f\t",gsl_matrix_get(gaussians->VV,dd1,dd1));
      for (dd1=0; dd1 != d-1; ++dd1)
	for (dd2=dd1+1; dd2 != d; ++dd2){
	  fprintf(logfile,"%f\t",gsl_matrix_get(gaussians->VV,dd1,dd2));
	}
      ++gaussians;
      fprintf(logfile,"\n\n");
    }
    gaussians -= K;
    fflush(logfile);
  }





  //Then run projected_gauss_mixtures
  proj_gauss_mixtures(data,N,gaussians,K,(bool *) fixamp,(bool *) fixmean, (bool *) fixcovar,avgloglikedata,tol,(long long int) maxiter, likeonly, w,splitnmerge);


  //Print the final model parameters to the logfile
  if (keeplog){
    fprintf(logfile,"\nFinal model parameters obtained:\n\n");
    for (kk=0; kk != K; ++kk){
      fprintf(logfile,"Gaussian ");
      fprintf(logfile,"%i",kk);
      fprintf(logfile,"\n");
      fprintf(logfile,"amp\t=\t");
      fprintf(logfile,"%f",(*gaussians).alpha);
      fprintf(logfile,"\n");
      fprintf(logfile,"mean\t=\t");
      for (dd1=0; dd1 != d; ++dd1){
	fprintf(logfile,"%f",gsl_vector_get(gaussians->mm,dd1));
	if (dd1 < d-1) fprintf(logfile,"\t");
      }
      fprintf(logfile,"\n");
      fprintf(logfile,"covar\t=\t");
      for (dd1=0; dd1 != d; ++dd1)
	fprintf(logfile,"%f\t",gsl_matrix_get(gaussians->VV,dd1,dd1));
      for (dd1=0; dd1 != d-1; ++dd1)
	for (dd2=dd1+1; dd2 != d; ++dd2){
	  fprintf(logfile,"%f\t",gsl_matrix_get(gaussians->VV,dd1,dd2));
	}
      ++gaussians;
      fprintf(logfile,"\n\n");
    }
    gaussians -= K;
    fflush(logfile);
  }



  //Then update the arrays given to us by IDL, but first, let's make sure wer're returning the right thing, ie the upper right of the VV matrix
  for (jj = 0; jj != K; ++jj){
    *(amp++) = gaussians->alpha;
    for (dd1 = 0; dd1 != d; ++dd1)
      *(xmean++) = gsl_vector_get(gaussians->mm,dd1);
    for (dd1 = 0; dd1 != d; ++dd1)
      for (dd2 = dd1+1; dd2 != d ; ++dd2)
	gsl_matrix_set(gaussians->VV,dd2,dd1,gsl_matrix_get(gaussians->VV,dd1,dd2));
    for (dd1 = 0; dd1 != d; ++dd1)
      for (dd2 = 0; dd2 != d; ++dd2)
	*(xcovar++) = gsl_matrix_get(gaussians->VV,dd1,dd2);
    ++gaussians;
  }
  gaussians -= K;
  amp -= K;
  xmean -= K*d;
  xcovar -= K*d*d;
  
  //And free any memory we allocated
  for (ii = 0; ii != N; ++ii){
    gsl_vector_free(data->ww);
    gsl_matrix_free(data->SS);
    gsl_matrix_free(data->RR);
    ++data;
  }
  data -= N;
  free(data);
  
  for (jj = 0; jj != K; ++jj){
    gsl_vector_free(gaussians->mm);
    gsl_matrix_free(gaussians->VV);
    ++gaussians;
  }
  gaussians -= K;
  free(gaussians);

  if (keeplog)
    fclose(logfile);

  return 3;
}





bool cleanup(){
  return true;
}

int main (void){

  //Initialize some of the arrays (will be reallocated once more info is available)


  //Loop over commandline arguments to extract the different filenames


  //Read initial conditions file
  printf ("Reading initial conditions file...\n");
  if ( read_IC("IC.dat") )
    printf ("Successfully read initial conditions file\n");
  else {
    printf ("Reading initial conditions failed\n");
    printf("Returning");
    return -1;
  }

  //Successfully implemented reading the initial conditions file

  //Next up: reading the datafile
  //put the data in a structure that has [d_i,ww_i,SS_i,PP_i]
  //[d_i is a byte array which contains which values are missing


  //Read datafile  
  printf("\nReading data...\n");
  if (read_data("inputdata.dat"))
    printf("Succesfully read the data\n");
  else {
    printf("Couldn't read the data...\n");
    printf("Returning\n");
    return -1;
  }


  //Try showing the data
  int ii;
  for (ii=0; ii != N-1; ++ii){
    printf("datapoint %i:\n",ii+1);
    printf("Mean\t\t\t=\t%f\t%f\n",gsl_vector_get(data->ww,0),gsl_vector_get(data->ww,1));
    fflush(stdout);
    printf("Variance\t\t=\t%f\t%f\t%f\n",gsl_matrix_get(data->SS,0,0),gsl_matrix_get(data->SS,1,1),gsl_matrix_get(data->SS,0,1));
    printf("Projection matrix\t=\t%f\t%f\t%f\t%f\n",gsl_matrix_get(data->RR,0,0),gsl_matrix_get(data->RR,0,1),gsl_matrix_get(data->RR,1,0),gsl_matrix_get(data->RR,1,1));
    ++data;
  }   
  for (ii=N-1; ii != N; ++ii){
    printf("datapoint %i:\n",ii+1);
    printf("Mean\t\t\t=\t%f\n",gsl_vector_get(data->ww,0));
    fflush(stdout);
    printf("Variance\t\t=\t%f\n",gsl_matrix_get(data->SS,0,0));
    printf("Projection matrix\t=\t%f\t%f\n",gsl_matrix_get(data->RR,0,0),gsl_matrix_get(data->RR,0,1));
    ++data;
  } 
  data -= N;
  //data = startdata;



  //Successfully implemented reading the data into a structure!!

  //Next up: ouputting the result to a file (should basically be the same as the code at the end of read_IC
  //Done!


  //Then decide on the function EM
  //Run EM
  double avgloglikedata;
  proj_gauss_mixtures(data, N, gaussians, K,fixampP, fixmeanP, fixcovarP, &avgloglikedata, tol, maxiter, likeonly, w,splitnmerge);





  //Write the result
  printf("Writing the results...\n");
  if (write_model("result.dat") )
    printf("Succesfully wrote the results...\n");
  else {
    printf("Couldn't write the results...\n");
    printf("Returning\n");
    return -1;
  }


  //clean up
  if (cleanup());
  else {
    printf("Couldn't quite clean everything up...\n");
    printf("Sorry!\n");
    return -1;
  }

  return 0;
}
