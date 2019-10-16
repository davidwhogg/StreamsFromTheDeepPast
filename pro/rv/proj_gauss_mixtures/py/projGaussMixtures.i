/* projGaussMixtures.i */
%module projGaussMixtures
%{
#define SWIG_FILE_WITH_INIT
extern int proj_gauss_mixtures_IDL(double * ydata, double * ycovar, 
                            double * projection, int N, int dy, 
                            double * amp, double * xmean, 
                            double * xcovar, int d, int K, 
                            char * fixamp, char * fixmean, 
                            char * fixcovar, 
                            double * avgloglikedata, double tol, 
                            int maxiter, char likeonly, double w, 
                            char * logfilename, int slen, int splitnmerge,
                            char * convlogfilename, int convloglen);
%}

%include "carrays.i"
%array_class(int, intArray);
%array_class(double, doubleArray);
%apply (char *LOGFILENAME, int SLEN) { (char *logfilename, int slen) };
%apply (char *CONVLOGFILENAME, int CONVLOGLEN) { (char *convlogfilename, int convloglen) };
%rename proj_gauss_mixtures_IDL proj_gauss_mixtures;
extern int proj_gauss_mixtures_IDL(double * ydata, double * ycovar, 
                            double * projection, int N, int dy, 
                            double * amp, double * xmean, 
                            double * xcovar, int d, int K, 
                            char * fixamp, char * fixmean, 
                            char * fixcovar, 
                            double * avgloglikedata, double tol, 
                            int maxiter, char likeonly, double w, 
                            char * logfilename, int slen, int splitnmerge,
                            char * convlogfilename, int convloglen);

