#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fftw3.h>
int main(int argc , char **argv)
{
     double *in;
     double *tmpin;
     fftw_complex *out;
     fftw_complex *tmpout;
     fftw_plan p;
     char buf[512];
     FILE*file=NULL;
     int k,lines=0;
     char *ret;
     int data,average=0,stride=0;
     double dt,val;


     if(argc==1) exit(printf("Usage: %s dt [maxlines average stride ]\n",argv[0]));
     file=stdin;
     dt = atof(argv[2]);
     if(argc>=3)  lines=atoi(argv[2]);
     if(argc>=4)  average = atoi(argv[3]); 
     if(argc==5)  stride = atoi(argv[4]);


     if(lines==0) { data=1024; } else { data=lines; } ; 
     in = (double*)fftw_malloc(data*sizeof(double));
     out= (fftw_complex*)fftw_malloc(data*sizeof(fftw_complex));
     k=0;
     while( (ret=fgets(buf,512,file))!=NULL){
        if(buf[0]!='#'){ 
	   sscanf(buf,"%lf",&val);
           in[k]=val;
           k++;
	   if(k>=data) {
		data*=2;
		tmpin=(double*)fftw_malloc(data*sizeof(double));
		tmpout=(fftw_complex*)fftw_malloc(data*sizeof(fftw_complex));
     		memcpy(tmpin, in, (k)*sizeof(double));
     		memcpy(tmpout,out, (k)*sizeof(fftw_complex));
		free(in);
		fftw_free(out);
		in=tmpin;
		out=tmpout;
	   }
        }
     }
     in[0]/=2; 
     data=k;
     tmpin=(double*)fftw_malloc(data*sizeof(double));
     tmpout=(fftw_complex*)fftw_malloc(data*sizeof(fftw_complex));
     memcpy(tmpin, in, (data)*sizeof(double));
     memcpy(tmpout,out, (data)*sizeof(fftw_complex));
     free(in);
     fftw_free(out);
     in=tmpin;
     out=tmpout;





     fprintf(stderr,"lines=%d :avg = %d stride = %d\n",data,average,stride);
     p = fftw_plan_dft_r2c_1d(data, in, out,FFTW_ESTIMATE);
     fftw_execute(p);
     if(argc>4){
       average--;
       while (average>0) { 
        fftw_complex * tmpout = (fftw_complex*)fftw_malloc((data-average*stride)*sizeof(fftw_complex));
        p = fftw_plan_dft_r2c_1d(data-average*stride, in, tmpout,FFTW_ESTIMATE);
	fftw_execute(p);
        fprintf(stderr,"%d : k=0 -> %f\n",average,tmpout[0][0]);
        for(k=0;k<data-average*stride;k++){
		out[k][0]+=tmpout[k][0];
		out[k][1]+=tmpout[k][1];
        }
	fftw_free(tmpout);
        average--;
       }
       average=atoi(argv[2]); 
     } else { average = 1 ; }
     for(k=0;k<data-average*stride;k++){
	printf("%f %g %g\n",k*2*M_PI/(dt*data),out[k][0]/average,out[k][1]/average);
     }
      
} 
//Programs using RFFTW should link with -lrfftw -lfftw -lm on Unix, or with the FFTW, RFFTW, and math libraries in general.

