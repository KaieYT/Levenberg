#include "stdio.h"
#include "stdlib.h"
#define NR_END 1
#define FREE_ARG char*

double *Pvector(long size)
{
	double *v;

	v = (double *)malloc((size_t)((size) * sizeof(double)));
	if(v == NULL)	fprintf(stderr, "Error : allocation failure in vector().\n");
	return v;
}

void free_vector(double *v)
{
	free(v);
}
void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
void free_ivector(int *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double **matrix(long rows, long cols)
{
	double **m, **p;
	long i;

	m = (double **)malloc((size_t)(rows * sizeof(double *)));
	if(m == NULL)
		printf("Error : allocation failure 1 in matrix().\n");
	*m = (double *)malloc((size_t)(rows * cols * sizeof(double)));
	if(*m == NULL)
		printf("Error : allocation failure 2 in matrix().\n");
	for(p = m, i = 1; i < rows; i++, p++)
		*(p + 1) = *p + cols;
	return m;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
    if (!m){
        printf("allocation failure 1 in matrix()");
        exit(1);
    }
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
    if (!m[nrl]){
        printf("allocation failure 2 in matrix()");
        exit(1);
    }
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float *Yvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if (!v){
        printf("allocation failure in vector()");
        exit(1);
    }
	return v-nl+NR_END;
}
int *ivector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
    if (!v){
        printf("allocation failure in vector()");
        exit(1);
    }
	return v-nl+NR_END;
}

void free_matrix(double **m)
{
	free(m);
}
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

