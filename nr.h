
double *Pvector(long size);
void free_vector(double *v);
double **matrix(long rows, long cols);
float **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
float *Yvector(long nl, long nh);
int *ivector(long nl, long nh);

void nrerror(char error_text[]);

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
/*
double f1dim(double x);
double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
void  powell(double *p, double **xi, int n, double ftol, int *iter, double *fret, double (*func)(double *));
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []));
int mainPowell();
*/

