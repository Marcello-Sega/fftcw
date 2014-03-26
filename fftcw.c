#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <fftw3.h>
#define DEBUG_LEVEL 0
typedef enum {MEM_ADD,MEM_RESIZE,MEM_REMOVE} mem_operation;
typedef enum {TYPE_REAL,TYPE_COMPLEX} number_type;
int memory_usage=0;
void usage (char*);
void memusage(int amount,char* address,mem_operation op);
void compute_corr_ft(int nobj, int dim, int window, int nframes, fftw_complex * out, double * autocorr, double * croscorr);
int suck (int nobj, int dim, int maxframes,  double **buffer,FILE *cid);
int do_fft(int nobj, int dim, int window, int nframes, int totframes, double ** buffer,double ** autocorr, double ** croscorr );
void do_spectrum(int nobj, int dim, int window, int nframes, int cut, double * in, fftw_complex ** out);
void do_back_fft(int nobj, int dim, int window, int nframes, int nfiles, double ** autocorr, double ** croscorr, double ** p_autocorr_direct, double ** p_croscorr_direct);
int output (int nobj, int window, int nframes, double dt, void * ac_dir, void * cc_dir, FILE*outfile, number_type type,int RSpaceOutputCut,int KSpaceOutputCut);

int main(int argc, char ** argv) { 

    double * buffer=NULL;
    double timestep=1;
    char ch;
    int dim=0,nobj=0,nframes=0,totframes,window=0,maxframes=0,cut=0;
    int i,j,k;
    double * autocorr, * autocorr_direct, * croscorr, * croscorr_direct;
    fftw_complex  * autocorr_spectrum=NULL, * croscorr_spectrum=NULL;
    char filename [1024];
    char awkstring [1024];
    FILE * out;
    int RSpaceOutputCut=-1,KSpaceOutputCut=-1;
    int nfiles=0;
    char infilenames[512][1024];
    FILE * in=NULL;

    extern char *optarg;
    extern int optind;

    if(getenv("EXT_COMM")) { 
	strcpy(awkstring,getenv("EXT_COMM"));
        printf("Using the external command: %s to pre-treat data\n ",awkstring);
    } else { 
	strcpy(awkstring,"-");
    }

    sprintf(filename,"-");
    i=j=k=0;
    while ((ch = getopt (argc, argv, "n:d:f:o:w:m:t:c:T:x:R:F:")) != -1)
    {
        switch (ch)
        {
             case 'T': 
                optind--;
                for(;optind<argc && argv[optind][0]!='-'; optind++){
			sprintf(infilenames[nfiles],"%s",argv[optind]); 
                        nfiles++;
                        if(nfiles==512) exit(printf("Error, hardcoded max number of input files (512) reached\n"));
                     }
                break;
	     case 'R': RSpaceOutputCut = atoi(optarg); i++; break;
	     case 'F': KSpaceOutputCut    = atoi(optarg); i++; break;
             case 'd': dim    = atoi(optarg); i++; break;
             case 'n': nobj   = atoi(optarg); i++; break;
             case 'w': fprintf(stderr,"warning: -w option unused yet\n");break; 
             case 'm': maxframes= atoi(optarg);    break; 
             case 'f': nframes= atoi(optarg);      break;
             case 'c': cut= atoi(optarg);      break;
             case 'o': strncpy(filename,optarg,1024); break;
             case 't': timestep=atof(optarg); i++; break;
             case 'h': case '?': default: usage (argv[0]);
        }
    }
    if(i<3) usage(argv[0]);
    if(nfiles==0) { nfiles=1; sprintf(infilenames[0],"-"); }  
    for(i=0;i<nfiles;i++){
         if(in!=NULL) {
              if(!strcmp(awkstring,"-")) {
                fclose(in); 
                fflush(in);
              } else {
                pclose(in); 
                fflush(in);
              }
              in=NULL;
          }
         if(!strcmp(infilenames[i],"-")){ in=stdin; }  else { 
                  if(!strcmp(awkstring,"-")) {
                       in=fopen(infilenames[i],"r");  
                  } else { 
                       char command[8192]={""};
                       strcat(command,awkstring); strcat(command," "); strcat(command,infilenames[i]);
                       fprintf(stderr,"reading through: %s\n",command);
		       in = popen(command,"r");
                  } 
                  if(in==NULL) exit(printf("Problems reading file %s . Exiting.\n",infilenames[i]));
         }
         totframes=suck(nobj,dim,maxframes, &buffer,in);  /* Obtain the full dataset */
         if(RSpaceOutputCut==-1)RSpaceOutputCut=totframes;
         if(KSpaceOutputCut==-1)KSpaceOutputCut=totframes;
         if(nframes==0) nframes=totframes;
         if(cut==0) cut=nframes;
         fprintf(stderr,"Data acquired (file %s, %d samples), will compute correlations of length %d\n",infilenames[i],totframes,nframes);


#if (DEBUG_LEVEL >= 1)
    if(totframes>0){
            printf("Sucked %d frames (%d objects of %d elements)\n",totframes,nobj,dim);
    }else {exit(printf("Return value = %d\n",totframes));}
#if (DEBUG_LEVEL >= 3)
    for(i=0;i<totframes;i++){
        for(j=0;j<nobj;j++){
            for(k=0;k<dim;k++) printf("%f ",buffer[i*nobj*dim+j*dim+k]);
            printf("\n");
        }
    }
#endif 
#endif 

        /* 
           to be implemented : repeat for different starting point 
           (should follow NR overlap-add method 
           [Num. Rec C, page 544, 2nd ed. Cambridge]) 
        */
    
        if(window==0) window = nframes;
    
        /* Don't get why the code does not work with odd number so far, so let's restrict to even ones */
        if(window%2) window--;
        if(nframes%2) nframes--;
        
        fprintf(stderr,"Performing fwd FFTs: nobj=%d dim=%d window=%d nframes=%d buffer=%p\n",nobj,dim,window,nframes,(void*)buffer);
        do_fft(nobj,dim,window,nframes,totframes,&buffer,&autocorr,&croscorr); 
        free(buffer); buffer=NULL;
    }
    
    fprintf(stderr,"Performing bck FFTs\n");
    do_back_fft(nobj, dim, window, nframes, nfiles, &autocorr, &croscorr, &autocorr_direct, &croscorr_direct);
    
    
    if(!strcmp("-",filename)) { out=stdout; } 
    else {  out = fopen(filename,"w");} 
    output (nobj, window, nframes, timestep, (void*)autocorr_direct, (void*)croscorr_direct,out,TYPE_REAL,RSpaceOutputCut,KSpaceOutputCut);
    if(strcmp("-",filename)) fprintf(stderr,"Time correlations saved in file %s\n",filename);
#if 0 
    fprintf(stderr,"Computing Spectra (without cut)\n");
    
    do_spectrum(nobj,dim,window,nframes,window,autocorr_direct,&autocorr_spectrum); 
    if(nobj>1){ /* then there are xcorrelations, otherwise,  no */
        do_spectrum((nobj*(nobj-1))/2,dim,window,nframes,window,croscorr_direct,&croscorr_spectrum); 
    }

    if(!strcmp("-",filename)) { out=stdout; } 
    else { sprintf(filename,"%s.fft",filename) ; out = fopen(filename,"w");} 
     /* do_spectrum() computes the fft from a reduced window of 'cut' elements, which we need to pass as
        a window value to output() [actually, cut/2+1 ]*/
    output (nobj, window/2+1 , nframes, timestep, (void*)autocorr_spectrum, (void*)croscorr_spectrum,out,TYPE_COMPLEX,RSpaceOutputCut,KSpaceOutputCut);
    if(strcmp("-",filename)) fprintf(stderr,"Spectra saved in file %s\n",filename);

    free(croscorr_spectrum); croscorr_spectrum=NULL;
    free(autocorr_spectrum); autocorr_spectrum=NULL;
    fprintf(stderr,"Computing Spectra (with cut on the autocorrelation)\n");
    do_spectrum(nobj,dim,window,nframes,cut,autocorr_direct,&autocorr_spectrum); 
    if(nobj>1){ /* then there are xcorrelations, otherwise,  no */
        do_spectrum((nobj*(nobj-1))/2,dim,window,nframes,cut,croscorr_direct,&croscorr_spectrum); 
    }

    if(!strcmp("-",filename)) { out=stdout; } 
    else { sprintf(filename,"%s.cut",filename) ; out = fopen(filename,"w");} 
     /* do_spectrum() computes the fft from a reduced window of 'cut' elements, which we need to pass as
        a window value to output() [actually, cut/2+1 ]*/
    output (nobj, cut/2+1 , nframes, timestep, (void*)autocorr_spectrum, (void*)croscorr_spectrum,out,TYPE_COMPLEX,RSpaceOutputCut,KSpaceOutputCut);
    if(strcmp("-",filename)) fprintf(stderr,"Spectra with cut saved in file %s\n",filename);
#endif

    
    return 1;
} 

void usage (char *arg){ exit(printf( "Usage : %s -t <timestep> -n <n_objects> -d <dimension>  [ -w <window> ] [ -f <nframes> ] [ -o <output> ] [ -m <maxframes> ] [ -c <spectrum_input_cut> ] [ -T <inputfile> [ <inputfile2> [...] ] ] [ -F <spectrum_output_cut> ] [ -R <correlation_output_cut> ] [ -h ]\nNote: you can supply a command through which to pipe your data before processing them using the EXT_COMM variable\n e.g. EXT_COMM=\"/usr/bin/awk \'{print \\$1}\'\"\n",arg  ));}

void memusage(int amount,char * address,mem_operation op){
#if DEBUG_LEVEL <=-1
     static int last=0;
     int increase;
     if(address==NULL) exit(printf("Problems allocating memory (@ %p)\n",address));
     switch(op){ 
          case MEM_ADD:    last=increase=amount; fprintf(stderr,"malloc: ");break;
          case MEM_RESIZE: increase=amount-last;last=amount; fprintf(stderr,"realloc: "); break;
          case MEM_REMOVE: last=0; increase=-amount;  fprintf(stderr,"free: "); break;
          default: return ;
     }
     memory_usage += increase;
     fprintf(stderr,"Memory usage %.1f kB (%.1f kB increase @ %p)\n",(double)memory_usage/8/1024,(double)increase/8/1024.,address);
#endif
}

int suck (int nobj, int dim, int maxframes, double **buffer, FILE * cid){

    static char * str=NULL, *ret;
    double *buf=NULL;
    char * token;
    int strsize; 
    int bufsize=1024,frames,j,k;
    strsize = 256*dim;

    if(str==NULL) { str = (char*)malloc(strsize*sizeof(char)); memusage(strsize*sizeof(char),str,MEM_ADD); }  
      if(str==NULL) exit(printf("Error allocating memory\n"));
    if(buf==NULL) buf=(double*) malloc(bufsize * dim * nobj * sizeof(double));
    memusage(bufsize * dim * nobj * sizeof(double),(char*)buf,MEM_ADD);
    frames=0; 
    while (1) { 
        if(frames==bufsize) {
            /* let's expand our buffer exponentially, we will trim it at the end */
            fprintf(stderr,"%d frames read\n",frames);
            bufsize*=2;
            buf=(double*)realloc(buf,bufsize * dim * nobj * sizeof(double));
            memusage(bufsize * dim * nobj * sizeof(double),(char*)buf,MEM_RESIZE);
        }
        for(j=0 ; j < nobj; j++){
            ret=fgets(str,strsize,cid);
            /* EOF reached */
            if(ret==NULL) { 
                    buf=(double*)realloc(buf,frames * dim * nobj * sizeof(double)) ; 
                    memusage(frames* dim * nobj * sizeof(double),(char*)buf,MEM_RESIZE);
                    *buffer=buf;
                    return frames; 
            }
            if(strlen(ret)==1){ 
                    frames-=frames%nobj ;  
                    buf=(double*)realloc(buf, frames  * dim * nobj* sizeof(double)) ; 
                    memusage(frames * dim * nobj * sizeof(double),(char*)buf,MEM_RESIZE);
                    *buffer=buf;
                    return frames; 
            }/* SAW: check this, in case of malformed input ...  */

            /* This is needed, because fgets() retains the '\n', which would be then reported as a token */
            str[strlen(ret)-1]='\0';
               token = strtok (str, " ");
            if(token==NULL){ *buffer=buf;  return -1; } 
#if DEBUG_LEVEL >= 3 
            printf("-> %f\n",atof(token)); fflush(stdout);
#endif
            buf[frames*nobj*dim + j*dim]= atof(token);
            for(k=1 ; k < dim ; k++) {
                token = strtok (NULL," ");
                buf[frames*nobj*dim + j*dim + k ] = atof(token);
            }
        }
        frames++;
        if(frames==maxframes) { 
             buf=(double*)realloc(buf,frames * dim * nobj * sizeof(double)) ; 
             memusage(frames * dim * nobj * sizeof(double),(char*)buf,MEM_RESIZE);
             *buffer=buf; 
             return frames;
        }
    }
}

void do_spectrum(int nobj, int dim, int window, int nframes, int cut, double * in, fftw_complex ** out)
{
   int i;
   fftw_plan pfwd=NULL;
   int n[1];
   int fft_size;
   int rank=1;
   int howmany;
   n[0]=cut;
   for(i=0;i<nobj;i++) in[i]/=2.; /* let's make the Laplace-Fourier Transform (let's put that back later on...)*/ 
   fft_size=((n[0]*nobj))*sizeof(fftw_complex);
   /* TODO: not reusable so far...handle reallocation */
   howmany=nobj;
   
   *out = (fftw_complex*) fftw_malloc(fft_size);
   memusage(fft_size,(char*)(*out),MEM_ADD);

   if(pfwd==NULL) pfwd = fftw_plan_many_dft_r2c(rank, n, howmany, in, n, howmany, 1 , *out, n, howmany , 1 , FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
   if (pfwd==NULL) exit(fprintf(stderr,"Something's wrong with planning fftw (line %d)\n",__LINE__));

   fftw_execute(pfwd);  
   memusage(nobj*(window+nframes)*sizeof(double),(char*)in,MEM_REMOVE);
   fftw_destroy_plan(pfwd); pfwd=NULL;
   for(i=0;i<nobj;i++) in[i]*=2.;
   return ;
}

int do_fft(int nobj, int dim, int window, int nframes, int totframes, double ** buffer,double ** p_autocorr, double ** p_croscorr){

    /* notation the same as in fftw3 manual, for clarity */
    int rank    = 1;    /* we are computing 1d transforms */
    int n[1];           /* data + pad */
    int howmany;        /* dim * nobj */
    int idist=1, odist=1;
    int istride, ostride;    /* distance between two elements in
                                       the same column */
    int *inembed , *onembed;
    static fftw_plan pfwd = NULL;
    static double * in = NULL;
    fftw_complex * out = NULL;
    static double * autocorr = NULL;
    static double * croscorr = NULL;

    int i;


    howmany     =    dim * nobj;
    n[0]        =    nframes+window; /* the second part to be zero-padded */
    istride     =    ostride = howmany;    
    inembed     =    onembed = n;
    
    if(in==NULL){
         in = (double*) fftw_malloc(sizeof(double) * howmany * n[0]);
         memusage(sizeof(double) * howmany * n[0],(char*)in,MEM_ADD);
    }

    /* this is for (to be implemented) windowing...TODO: redesign */ 
    memcpy(in,*buffer,sizeof(double) * howmany * nframes);
    for(i = howmany * nframes ; i <  howmany * n[0]  ; i++ ) in[i]=0.0; /* zero padding */
    if(totframes==nframes){ /* then we can safely discard buffer and save some memory */
       free(*buffer); *buffer=NULL;
       memusage(sizeof(double) * howmany * nframes,(char*)*buffer,MEM_REMOVE);
    }
    if(out==NULL){ out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * howmany * n[0]); memusage(sizeof(fftw_complex) * howmany * n[0],(char*)out,MEM_ADD); }
    if(out==NULL) exit(printf("Error allocating memory\n"));
    if(autocorr==NULL){ 
                 autocorr = (double*) fftw_malloc(sizeof(double) * nobj * (n[0]/2+1)); 
                 memusage(sizeof(double) * nobj * (n[0]/2+1),(char*)autocorr,MEM_ADD); 
                 for(i=0;i< nobj * (n[0]/2+1);i++) autocorr[i]=0.;
    }
    if(croscorr==NULL){ 
                 croscorr = (double*) fftw_malloc(sizeof(double) * ((nobj*(nobj-1))/2) * (n[0]/2+1)); 
                 memusage(sizeof(double) *  ((nobj*(nobj-1))/2)  * (n[0]/2+1),(char*)autocorr,MEM_ADD); 
                 for(i=0;i< ((nobj*(nobj-1))/2) * (n[0]/2+1);i++) croscorr[i]=0.;
    }
    *p_autocorr = autocorr;
    *p_croscorr = croscorr;

#if (DEBUG_LEVEL >= 3)
    {
        int i,j,k;
        printf("actual data going to be fed to fftw:\n");
        for(i=0;i<n[0];i++){
            for(j=0;j<nobj;j++){
                for(k=0;k<dim;k++) printf("%f ",in[i*nobj*dim+j*dim+k]);
                printf("\n");
            }
        }
    }
#endif 
    /* we need to use FFTW_ESTIMATE, as we are computing only one (or few) FFTs! */
    for(i=0;i<n[0]*howmany;i++) { out[i][0]=out[i][1]=1.0 ; } 
    if(pfwd==NULL){
         pfwd = fftw_plan_many_dft_r2c(rank, n, howmany, in, inembed, istride, idist, out, onembed, 
                                     ostride, odist, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
            if (pfwd==NULL) exit(fprintf(stderr,"Something's wrong with planning fftw (line %d)\n",__LINE__));

    }

    fftw_execute(pfwd);    

    compute_corr_ft(nobj, dim,  window, nframes, out, autocorr, croscorr);

    fftw_free(out); 
    fftw_destroy_plan(pfwd); pfwd=NULL;
    memusage(n[0]*dim*nobj*sizeof(fftw_complex),(char*)out,MEM_REMOVE);
    out=NULL;
    return 0;
}

void compute_corr_ft(int nobj, int dim, int window, int nframes,  fftw_complex * out, double * autocorr, double * croscorr){
    int i,j,jj,k,index,index2,xobj;
   
    int n[1],nxcorr;
    n[0]=(nframes+window)/2+1;
    /* let's start by filling the FT'ed autocorrelations (real) ... */
    for    (i = 0 ; i < n[0]; i++){ 
      for(j = 0 ; j < nobj ; j++){ 
        index = i * nobj * dim  + j * dim ; 
        for(k = 0 ; k < dim ; k++){  /* we increment: anyway, the back fourier transform to be taken later is linear ... */
          autocorr[i * nobj + j] += (out[index + k][0] * out[index + k][0] + out[index + k][1] * out[index + k][1])/(nframes+window);
#if DEBUG_LEVEL >=2 
	  printf("row=%d obj=%d dir=%d FT= %f , %f  Conv = %f\n", i,j,k,out[index + k][0], out[index + k][1],autocorr[i * nobj + j]);
#endif
        }
      }    
    }
    /* ... and continue with the FT'ed cross correlations */

    if(nobj-1==0) return; /* there are no cross correlations to be calculated*/
    /* NOTE: we take 1/2 [ corr(a , b) + corr(b , a) ] (which is equivalent to use time-reversal symmetric signal)  
             and which leads to computing F(a) conj(F(b)) + F(b) conj(F(a)). This has no imaginary part, and we can use real DFTs... */ 
    nxcorr=(nobj*(nobj-1))/2;
    for (i = 0 ; i < n[0]; i++){ 
      xobj=0;
      for(j = 0 ; j < nobj ; j++){ 
        index = i * nobj * dim  + j * dim ; 
        for(jj = j+1 ; jj < nobj ; jj++){ 
            index2 = i * nobj * dim  + jj * dim ; 
            for(k = 0 ; k < dim ; k++){  /* we increment: anyway, the back fourier transform to be taken later is linear ... */
                croscorr[i * nxcorr + xobj ]+= (out[index + k][0] * out[index2 + k][0] + out[index + k][1] * out[index2 + k][1])/(nframes+window);
#if DEBUG_LEVEL >=2 
	    printf("row=%d obj1=%d obj2=%d xobj=%d dir=%d FTs=(%f * %f + %f * %f)  XConv = %f\n", i,j,jj,xobj,k,
			out[index + k][0] , out[index2 + k][0] , out[index + k][1] , out[index2 + k][1],
			croscorr[i * nxcorr+ xobj]);
#endif
            }
            xobj++;
        }
      }    
    }
}


void do_back_fft(int nobj, int dim, int window, int nframes, int nfiles, double ** autocorr, double ** croscorr, double ** p_autocorr_direct, double ** p_croscorr_direct)
{
/* This gives the autocorrelation function. A subsequent Laplace-FT will give the spectrum. */

/*  Note: by computing fa=fft(a) [with a zero-padded] ; corr=ifft(fa.*conj(fa))
 *  corr is so that every element has a statistics of size(a). This is however
 *  is induced by the padding, as correlations of lag 0 have statistics of size(a),
 *  but correlations of lag 1 have statistics of size(a)-1, and so on. Hence, the 
 *  need to normalize properly [ 1./(nframes-i)  factor ] below.
 */ 


    static fftw_plan pbckAuto = NULL;
    static fftw_plan pbckCros = NULL;

    static double * autocorr_direct=NULL;
    static double * croscorr_direct=NULL;
    int i,j;
    /* notation the same as in fftw3 manual, for clarity */
    int rank    = 1;    /* we are computing 1d transforms */
    int n[1];           /* data size, no new padding to be added */
    fftw_r2r_kind kind[1];
    int howmany;        
    int nxcorr;
    int idist=1, odist=1;
    int istride, ostride;    /* distance between two elements in
                                       the same column */
    int *inembed , *onembed;

    /* Autocorrelations first */
    howmany = nobj ;
    n[0] = (nframes+window)/2+1;

    for(i=0;i<n[0]*howmany;i++) (*autocorr)[i] /= nfiles;

    istride     =    ostride = howmany;    
    inembed     =    onembed = n;
    kind[0]= FFTW_REDFT00;
    if (autocorr_direct==NULL){
           autocorr_direct = fftw_malloc(sizeof(double) * n[0] * howmany);
           memusage(sizeof(double) * n[0] * howmany,(char*)autocorr_direct,MEM_ADD);
    }
    * p_autocorr_direct = autocorr_direct;

    if (pbckAuto==NULL) pbckAuto = fftw_plan_many_r2r(rank, n, howmany, *autocorr, inembed, istride, idist, autocorr_direct, 
                                                      onembed, ostride, odist, kind, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
    if (pbckAuto==NULL) exit(fprintf(stderr,"Something's wrong with planning fftw (line %d)\n",__LINE__));

    fftw_execute(pbckAuto);
    fftw_free(*autocorr);
    memusage(sizeof(double) * n[0] * howmany,(char*)*autocorr,MEM_REMOVE);
    *autocorr=NULL;

    /* normalization */
    for(j=0;j<nobj;j++){
      for(i=0;i<window;i++){
		autocorr_direct[i*nobj+j]/=(nframes-i);
      }
    }

    
    /* and then cross correlations */
    howmany = (nobj*(nobj-1))/2;
    if(howmany==0) return; /* there are no cross correlations to be calculated*/
    n[0] = (nframes+window)/2+1;

    for(i=0;i<n[0]*howmany;i++) (*croscorr)[i] /= nfiles;

    rank=1;
    idist=odist=1;
    istride     =    ostride = howmany;    
    inembed     =    onembed = n;
    kind[0]= FFTW_REDFT00;
    if (croscorr_direct==NULL) { 
          croscorr_direct = fftw_malloc(sizeof(double) * n[0] * howmany);
          memusage(sizeof(double) * n[0] * howmany,(char*)croscorr_direct,MEM_ADD);
    }
    * p_croscorr_direct = croscorr_direct;

    if (pbckCros==NULL) pbckCros = fftw_plan_many_r2r(rank, n, howmany, *croscorr, inembed, istride, idist, croscorr_direct, 
                           onembed, ostride, odist, kind, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
    if (pbckCros==NULL) exit(fprintf(stderr,"Something's wrong with planning fftw (line %d)\n",__LINE__));
    fftw_execute(pbckCros);
    
    fftw_free(*croscorr); 
    memusage(sizeof(double) * n[0] * howmany,(char*)*croscorr,MEM_REMOVE);
    *croscorr=NULL;

    nxcorr=(nobj*(nobj-1))/2;
    
    /* normalization */
    for(j=0;j<nxcorr;j++){
      for(i=0;i<window;i++){
		croscorr_direct[i*nxcorr+j]/=(nframes-i);
      }
    }

    return;

}

int output (int nobj, int window, int nframes, double dt, void * ac_dir, void * cc_dir, FILE*outfile, number_type type, int RSpaceOutputCut, int KSpaceOutputCut){
       
    int i,j,jj,ind,nxcorr;
    double * pr_ac=NULL, * pr_cc=NULL;
    fftw_complex * pc_ac=NULL, * pc_cc=NULL;
    
    switch(type) {
	case TYPE_REAL: 
                 pr_ac=(double*) ac_dir; 
                 pr_cc=(double*) cc_dir; 
		 break;
        case TYPE_COMPLEX:
                 pc_ac=(fftw_complex*) ac_dir; 
                 pc_cc=(fftw_complex*) cc_dir; 
		 break;
        default: exit(printf("Type not implemented\n"));
    } 
    nxcorr=(nobj*(nobj-1))/2;
    for(j=0;j<nobj;j++){
      for(i=0;i<window;i++){
           if(type==TYPE_COMPLEX) {
		  if(i>=KSpaceOutputCut) break;
                  fprintf(outfile, "%d %d %f %f %f\n",j,j,2*M_PI*i/(window * dt),pc_ac[i*nobj+j][0],pc_ac[i*nobj+j][1]);
           } else {
		  if(i>=RSpaceOutputCut) break;
                  fprintf(outfile, "%d %d %f %f\n",j,j,i*dt,pr_ac[i*nobj+j]);
           }
      }
      fprintf(outfile,"\n");
    }
    ind=0;
    if(nxcorr>=1){
       for(j=0;j<nobj;j++){
         for(jj=j+1;jj<nobj;jj++){
            for(i=0;i<window;i++){
              if(type==TYPE_COMPLEX) {
                 fprintf(outfile, "%d %d %f %f %f\n",j,jj,2*M_PI*i/(window * dt),pc_cc[ i * nxcorr + ind][0],pc_cc[ i * nxcorr + ind][1]);
              } else { 
                 fprintf(outfile, "%d %d %f %f\n",j,jj,i*dt,pr_cc[ i * nxcorr + ind]);
              }
            }
            fprintf(outfile,"\n");
            ind++;
         }

       }
    }

    return 1;
}
