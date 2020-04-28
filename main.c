#include	<stdlib.h>
#include	<stdio.h>
#include	<errno.h>
#include	<time.h>
#include	<ctype.h>
#include	<math.h>
#include	<string.h>
#include	<memory.h>
//#include	"malloc.h"

#define TWOPI	6.28318530717959
#define	MAX(A,B) (((A)>(B))?(A):(B))
#define	MIN(A,B) (((A)<(B))?(A):(B))
#define	EPS 0.000000001 /* A tiny fraction to help with roundoff errors */
#define	ASIZE 20000000 /* allowable number of source data array points */
#define	IGSIZE 4096 /* allowable number of interpolation grid array points */
#define	LIMDEV 0.03 /* % allowable mean step deviation without interpolation */

/* dimension arrays explicitly, since we will not know how much memory
might be needed in each array */
double	xjbs[ASIZE+1];
double	yjbs[ASIZE+1];
volatile int columnSelect = 2;

void nrerror(char *error_text);
double *dvector(int nl,int nh);
void free_dvector(double *v,int nl,int nh);
void spline(double *x,double *y,int n,double yp1,double ypn,double *y2);
void splint(double *xa,double *ya,double *y2a,int n,double x,double *y);
int readxy(char *string,double *x,double *y,int *counter);
void output(int harm, double fund, double real, double imag, double span);
void reverse_array(double *arr, int n);
void move_time(double *arr, int n);
void remove_head(double *arr, int narr, int nremov);
void remove_tail(double *arr, int narr, int nremov);

void nrerror(error_text)
char error_text[];
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
	}

double *dvector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
	}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
	free((char*) (v+nl));
	nh;
	}

void spline(x,y,n,yp1,ypn,y2)
double x[],y[],yp1,ypn,y2[];
int n;
{
	int i,k;
	double p,qn,sig,un,*u,*dvector();
	void free_dvector();

	u=dvector(1,n-1);
	if (yp1 > 0.99e30)
	y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
		}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
		}
	if (ypn > 0.99e30)
	qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
		}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
	y2[k]=y2[k]*y2[k+1]+u[k];
	free_dvector(u,1,n-1);
	}

void splint(xa,ya,y2a,n,x,y)
double xa[],ya[],y2a[],x,*y;
int n;
{
	int klo,khi,k;
	double h,b,a;
	void nrerror();

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
		}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	}

int readxy(string,x,y,counter)
char *string;
double x[], y[];
int *counter;
{
	char* line;
	int	data=0;

	char linestyle[256];

    switch ( columnSelect )	{
        case 1:
        case 2:
            sprintf(linestyle, "%%le %%le");
            break;
        case 3:
            sprintf(linestyle, "%%le %%*le %%le");
            break;
        case 4:
            sprintf(linestyle, "%%le %%*le %%*le %%le");
            break;
        case 5:
            sprintf(linestyle, "%%le %%*le %%*le %%*le %%le");
            break;
        default:
            fprintf(stderr,"Column Selection greater than 5 not supported\n");
            exit(1);
    }

	/* ignore leading space character if it exists */
	if ( !strncmp(string," ",1) ) { line = string+1; }
	else { line = string; }
	/* and now look for a character that signals a number */
	if (	!strncmp(line,"1",1) ||
		!strncmp(line,"2",1) ||
		!strncmp(line,"3",1) ||
		!strncmp(line,"4",1) ||
		!strncmp(line,"5",1) ||
		!strncmp(line,"6",1) ||
		!strncmp(line,"7",1) ||
		!strncmp(line,"8",1) ||
		!strncmp(line,"9",1) ||
		!strncmp(line,"0",1) ||
		!strncmp(line,"+",1) ||
		!strncmp(line,"-",1) ||
		!strncmp(line,".",1) )	{  /* data line */
//		if (2!=sscanf(line, "%le %le", &x[*counter], &y[*counter])) {
		if (2!=sscanf(line, linestyle, &x[*counter], &y[*counter])) {
			/* did not get two floats! */
			fprintf(stderr,"dft: error reading line %d (did not get 2 numbers after \"numeric\" start):\n%s", *counter, line);
			exit(1);
			}
		if (*counter && x[*counter]<=x[*counter-1]) { /* check monotonicity */
			fprintf(stderr,"dft: Non-monotonic time value!\n");
			exit(1);
			}
		data=1;
		(*counter)++; /* count the number of data lines */
		}
	else {
		printf("%s", line); /* not a data line... echo */
		}
	return data;
	}

void output(harm, fund, real, imag, span)
int harm;
double fund,real,imag,span;
{
	double mag, pha;

	/* scale so that 1Vpk wave ideally gives magnitude 1 */
	mag = sqrt((real*real+imag*imag)*8.0)/span;
	/* now we scale for 1VRMS to give magnitude 1 */
	/* delete next line to give pk scaling again  */
	mag /= sqrt(2.0);
	pha = 360.0*atan2(imag,real)/TWOPI;
	printf("%le %le %0.2lf\n", harm*fund, mag, pha);
	fprintf(stderr,"harmonic %d = %lf < %.2lf                \r", harm, mag, pha);
	}

void reverse_array(arr, n)
double arr[];
int n;
{
    double temp;
    int a = 0;
    int b = n-1;

    while(a < b)
    {
        temp = arr[a];
        arr[a] = arr[b];
        arr[b] = temp;
        a++;
        b--;
    }
}

void move_time(arr, n)
double arr[];
int n;
{
    for(int i = n-1; i > 0; i--)
    {
        arr[i] = fabs(arr[i]-arr[0]);

    }
    arr[0] = 0; // set to 0 to avoid rounding issues
}

void remove_head(arr, narr, nremov)
double arr[];
int narr, nremov;
{
    for(int i = nremov-1; i < narr; i++)
        arr[i-nremov] = arr[i];

    remove_tail(arr, narr, nremov);
}

void remove_tail(arr, narr, nremov)
double arr[];
int narr, nremov;
{
    for(int i = (narr-nremov); i < narr; i++)
        arr[i] = 0;
}


int main(argc,argv)	/*pgm to take DFT of potentially ungridded data */
int 	argc;
char 	*argv[];
{
	FILE	*fp, *fopen();
	long	nlines=0L;
	int	ndatl=0;
	char	cinline[256];
	int	strncmp();
	int	i, j, ncycles, newn;
	int	harmonics, hstart=0;
	char	optn;
	double 	meandev, nstep;
	double	span, wspan, fundamental;
	double	hann;
	double	real, imag, ang;
	double	lsin, lcos;
	double	thissin, thiscos;
	double	width, kay, slope, dtemp, smallest;
	unsigned int	nbytes;
	void*	memsiz;
	double*	splines;
	double*	ny;
	volatile int argI , argD, argU, argL, argH, argT, argF;
	argI = argD = argU = argL = argH = argT = argF = 0;
	int endIdx = 0;

	if ( argc < 4 )	{ /* ?? */
		fprintf(stderr,"DFTp for windows V1.06 20M Points, VGF Sep 2019, JBS Dec 2013, after Apr-1995 & Dec 2012 \n");
		fprintf(stderr,"%d parameters is illegal.\n", argc-1);
		fprintf(stderr,"Usage: dftp ipfile w0 #harmonics [#starting_harmonic [I|D|U|L|C n||F|H n|T n]] [>opfile]\n");
		fprintf(stderr,"This program reads in ASCII files with one time/signal pair of values per\n");
		fprintf(stderr,"line and takes the DFT of signal*(raised_cosine) at w0 and #harmonics of w0,\n");
		fprintf(stderr,"starting at the #starting_harmonic (i.e., skipping over #starting_harmonics).\n");
		fprintf(stderr,"The data need NOT be equispaced and need NOT be a power of 2, but\n");
		fprintf(stderr,"must be in (increasing) monotonic time order.\n");
		fprintf(stderr,"The I option forces interpolation onto a regular grid, while\n");
		fprintf(stderr,"the D option forces a DFT without interpolation.\n");
		fprintf(stderr,"The U option inhibits windowing.\n");
		fprintf(stderr,"The L option subtracts a ramp from the data whose slope and offset is\n");
		fprintf(stderr,"found from a linear regression of the data in the window.\n");
        fprintf(stderr,"The C option selects the data column specified by the trailing number\n");
        fprintf(stderr,"There must be a space between C and n, n can be one of [2, 3, 4, 5] \n");
        fprintf(stderr,"The F option flips the data reversing the time order of the data\n");
        fprintf(stderr,"the last datum is now at time 0 and the first is at the of the original last point\n");
        fprintf(stderr,"this is achieved by subtracting the final time value from all time values and taking the absolute\n");
        fprintf(stderr,"OPTION F WILL REVERSE THE DATA BEFORE OPTION H OR OPTION T REMOVE CYCLES\n");
        fprintf(stderr,"The H option removes n leading cycles from the data\n");
        fprintf(stderr,"There must be a space between H and n, n must be a positive integer\n");
        fprintf(stderr,"The T option removes n trailing cycles from the data\n");
        fprintf(stderr,"There must be a space between T and n, n must be a positive integer\n");
        fprintf(stderr,"Lines which do not start in the first column with a character that starts a\n");
		fprintf(stderr,"legal number (+, -, ., or a digit 1-9 and 0) are treated as comment lines and\n");
		fprintf(stderr,"copied to the output verbatim; otherwise the line is treated as a data line.\n");
		fprintf(stderr,"One leading space at the start of a line is tolerated, but data lines which\n");
		fprintf(stderr,"do not yield two legal white-space-separated floats generate a fatal error.\n");
		fprintf(stderr,"Output is frequency/magnitude/phase triples.\n");
		exit(1);
		}
	if ((fp = fopen(argv[1], "rb")) == NULL)	{
		fprintf(stderr,"Cannot open %s as input.\n", argv[1]);
		exit(1);
		}

	sscanf(argv[2],"%le", &fundamental);
	sscanf(argv[3],"%d", &harmonics);
	if (argc>=5) {
		sscanf(argv[4],"%d", &hstart);
		fprintf(stderr,"Computing harmonics %d-%d with fundamental of %lf\n", hstart, hstart+harmonics, fundamental);
		} else {
		fprintf(stderr,"Computing first %d harmonics with fundamental of %lf\n", harmonics, fundamental);
		}

	if (argc>=6) { /* there is an option character */
        for(int i = 5; i < argc; i++)
        {
            sscanf(argv[i],"%c", &optn);
            switch ( optn )	{
                case 'i':
                case 'I': /* force interpolation */
                    argI = 1;
                    break;
                case 'd':
                case 'D': /* force NO interpolation (Direct) */
                    argD = 1;
                    break;
                case 'u':
                case 'U': /* inhibit windowing */
                    argU = 1;
                    break;
                case 'l':
                case 'L': /* subtract ramp*/
                    argL = 1;
                    break;
                case 'c':
                case 'C': /* select column*/
                    columnSelect = atoi(argv[i+1]);
                    fprintf(stderr, "Column Selected: %i\n", columnSelect);
                    if(columnSelect == 1)
                        fprintf(stderr, "Warning: first column must be time, reading data from column 2 instead\n");
                    i++;
                    break;
                case 'f':
                case 'F': /* flip data*/
                    argF = 1;
                    fprintf(stderr, "Flipping data order\n");
                    break;
                case 'h':
                case 'H': /* remove leading cycle*/
                    argH = atoi(argv[i+1]);
                    fprintf(stderr, "Leading cycles to remove: %i\n", argH);
                    i++;
                    break;
                case 't':
                case 'T': /* remove end cycle*/
                    argT = atoi(argv[i+1]);
                    fprintf(stderr, "Trailing cycles to remove: %i\n", argT);
                    i++;
                    break;
                default:
                fprintf(stderr,"Illegal option `%c'!\n", optn);
                exit(1);
			}
		}
    } else {
    optn = 'X';
    }

	/* Readin procedure */
	while ( fgets(cinline,254,fp) != NULL ) {
		nlines++;	/* count lines */
		readxy(cinline, xjbs, yjbs, &ndatl);
		if (ndatl>ASIZE) {
			fprintf(stderr,"Array overflow (Max %d points, got %d).\n", ASIZE, ndatl);
			exit(1);
			}
		}
	/* #chars will differ on DOS because newline=2 chars to DOS, 1 to unix */
	fprintf(stderr,"Read in %ld lines, %d datalines...       \n", nlines, ndatl);

	span = xjbs[ndatl-1]-xjbs[0];
	/* see what size window to multiply across the data */
	ncycles = (int)(fundamental*span*(1.0+EPS)); /* i = # whole periods in span */
	wspan = ncycles*(1/fundamental); /* largest span which is a multiple of the period */
	fprintf(stderr,"Natural span is %lf; window span (%d cycles fit) is %lf\n", span, ncycles, wspan);
	if ( ncycles<=0 ) {
		fprintf(stderr,"Natural span too small... aborting.\n");
		exit(1);
		}
	if (wspan<(span/2.0)) {
		fprintf(stderr,"Can't happen error 1 !\n");
		exit(1);
		}


    if((argT + argH) >= ncycles) {
		fprintf(stderr,"Unable to remove %d cycles from %d total cycles\n", (argT + argH), ncycles);
		exit(1);
        }

    if (argF)
    {
        fprintf(stderr,"Reversing data order\n");
        reverse_array(xjbs, ndatl);
        reverse_array(yjbs, ndatl);
        move_time(xjbs, ndatl);
    }
    if(argT) {
        fprintf(stderr,"Removing %d cycles from tail\n", argT);
        double time_to_remove = ((double)argT * (1.0/fundamental)) - EPS;
        if(time_to_remove >= xjbs[ndatl-1])
            fprintf(stderr,"Error unable to remove that many cycles from tail\n");

        int nremov = 0;
        for(int i = ndatl -1; i > 0; i--) {
            if((xjbs[ndatl-1] - xjbs[i]) >= time_to_remove)
                break;

            nremov++;
        }
        nremov--; // leave an extra datapoint behind

        remove_tail(xjbs, ndatl, nremov);
        ndatl -= nremov;
    }
    if(argH)
    {
        fprintf(stderr,"Removing %d cycles from head\n", argH);
        double time_to_remove = ((double)argH * (1.0/fundamental)) - EPS;
        if(time_to_remove >= xjbs[ndatl-1])
            fprintf(stderr,"Error unable to remove that many cycles from head\n");
        int nremov = 0;
        for(int i = 0; i < ndatl; i++)
        {
            if((xjbs[i] - xjbs[0]) >= time_to_remove)
                break;

            nremov++;
        }
        nremov--; // leave an extra data point behind

        remove_head(xjbs, ndatl, nremov);
        ndatl -= nremov;
    }

    if(argH || argT) {
        fprintf(stderr,"Recalculating span after removing %d cycles\n", (argH + argT));
        span = xjbs[ndatl-1]-xjbs[0];
        /* see what size window to multiply across the data */
        ncycles = (int)(fundamental*span*(1.0+EPS)); /* i = # whole periods in span */
        wspan = ncycles*(1/fundamental); /* largest span which is a multiple of the period */
        fprintf(stderr,"New span is %lf; window span (%d cycles fit) is %lf\n", span, ncycles, wspan);
        if ( ncycles<=0 ) {
            fprintf(stderr,"New span too small... aborting.\n");
            exit(1);
            }
        if (wspan<(span/2.0)) {
            fprintf(stderr,"Can't happen error 1 !\n");
            exit(1);
            }
    }

	/* are the data nearly equispaced? */
	dtemp = span/(ndatl-1);
	for(i=1, meandev=0.0; i<ndatl; i++) {
		meandev += fabs(dtemp-(xjbs[i]-xjbs[i-1]));
		}
	meandev /= span;
	fprintf(stderr,"Irregularity = %lf... ", meandev);

	if ( (meandev>LIMDEV && !argD) || argI) {
		/* interpolate to an equispaced grid */
		for(i=1, smallest=wspan; i<ndatl; i++) {
			smallest = MIN(smallest,xjbs[i]-xjbs[i-1]);
			}
		dtemp = MAX(2*ncycles*harmonics,wspan/smallest);
		for(newn=2; newn<dtemp && newn<=IGSIZE; newn*=2) {;}
		fprintf(stderr,"interpolating to grid of %d points. \n", newn);

		/* spline interpolation onto new grid */
		nbytes = (ndatl+1)*sizeof(double);
                if (nbytes < (ndatl+1)*sizeof(double)) {
                    fprintf(stderr,"IGSIZE too big for malloc!\n");
                    fprintf(stderr,"(Silly hashdefs in source file)\n");
                    exit(1);
                }
		if ( (memsiz=malloc(nbytes))==NULL ) {
			fprintf(stderr,"Insufficient memory (splines)! (wanted %ld)\n", (long)nbytes);
			exit(1);
			};
		splines=(double*)memsiz;
		nbytes = (newn)*sizeof(double);
		if ( (memsiz=malloc(nbytes))==NULL ) {
			fprintf(stderr,"Insufficient memory (ny)! (wanted %ld)\n", (long)nbytes);
			exit(1);
			}
		ny=(double*)memsiz;
		spline(xjbs-1,yjbs-1,ndatl,1.1e30,1.1e30,splines);
		/* now fill new array with interpolated values */
		nstep = wspan/(newn-1.0);
		for(i=0; i<newn; i++) {
			splint(xjbs-1,yjbs-1,splines,ndatl,xjbs[0]+i*nstep,&ny[i]);
			}
		/* plug new numbers into old arrays */
		for(i=0; i<newn; i++) {
			xjbs[i] = xjbs[0]+i*nstep;
			yjbs[i] = ny[i];
			}
		ndatl=newn;
		}
	if ((meandev<LIMDEV && !argI) || argD) {
		fprintf(stderr,"skipping interpolation.\n");
		/* it may be that there is not a point exactly at xjbs[wspan]... */
		for(i=0; i<ndatl && xjbs[i]-xjbs[0]<wspan-EPS; i++) {;} /* find first point outside wspan */
		if (i==ndatl) {
			fprintf(stderr,"Cannot happen! (xjbs[%d]-xjbs[0]=%le)-(wspan)=%le\n", i-1, xjbs[i-1]-xjbs[0], xjbs[i-1]-xjbs[0]-wspan);
			exit(1);
        }

		xjbs[i] = xjbs[0]+wspan; /* put it on the edge exactly */
//		yjbs[i] = yjbs[0]);	/* and copy the first datum around */
		yjbs[i] = yjbs[i-1] + (xjbs[i]-xjbs[i-1])*(yjbs[i+1]-yjbs[i-1])/(xjbs[i+1]-xjbs[i-1]);	/* and linearly interpolate the ,issing yvalue*/
        endIdx = i;
    }

    if (argL) {
		fprintf(stderr,"subtracting ramp.\n");
        if(!endIdx)
        {
            for(i=0; i<ndatl && xjbs[i]-xjbs[0]<wspan-EPS; i++) {;} /* find first point outside wspan */
            if (i==ndatl) {
                fprintf(stderr,"Cannot happen! (xjbs[%d]-xjbs[0]=%le)-(wspan)=%le\n", i-1, xjbs[i-1]-xjbs[0], xjbs[i-1]-xjbs[0]-wspan);
                exit(1); }
            endIdx = i;
        }

        /* fit a linear regression y = mx + b*/
        /* https://www.statisticshowto.datasciencecentral.com/probability-and-statistics/regression-analysis/find-a-linear-regression-equation/ */
        double sx = 0;  // sum x
        double sxy = 0; // sum xy
        double sxs = 0; // sum x^2
        double sy = 0;  // sum y
        double n = 0;   // number points

        for(int i = 0; i <= endIdx; i++)
        {
            sx += xjbs[i];
            sxy += (xjbs[i]*yjbs[i]);
            sxs += (xjbs[i]*xjbs[i]);
            sy += (yjbs[i]);
            n += 1;
        }

        double m = (n*sxy - sx*sy)/(n*sxs - sx*sx);
        double b = (sy*sxs- sx*sxy)/(n*sxs-sx*sx);
        fprintf(stderr, "\nM: %.8f B: %.8f\n", m, b);
        for(int i = 0; i <= endIdx; i++)
        {
            yjbs[i] = yjbs[i] - (m*xjbs[i] + b);	/* subtract the calculated slope from the data in the window*/
        }

    }

	/* apply a raised cosine window */
	if (!argU) {
		fprintf(stderr,"Windowing data...       \r");
		for(i=0; i<ndatl; i++) {
			if (xjbs[i]-xjbs[0] >= wspan) { hann=0.0; }
			else {	ang = TWOPI*(xjbs[i]-xjbs[0])/wspan;
				hann = 1.0-cos(ang);
				}
			yjbs[i] *= hann;
			}
		} else {
		fprintf(stderr,"Windowing inhibited---computing spectrum of unwindowed data.\n");
		}

	/* compute the DFTs, one for each requested harmonic */
	fprintf(stderr,"Transforming data...       \r");
	for(j=1+hstart; j<=harmonics+hstart; j++) {
		real = 0.0;
		imag = 0.0;
		lsin = sin(0.0);
		lcos = cos(0.0);
		kay = TWOPI * j * fundamental;
		for(i=1; i<ndatl; i++) {
			ang = (xjbs[i]-xjbs[0])*kay;
			width = xjbs[i]-xjbs[i-1];
			slope = (yjbs[i]-yjbs[i-1])/width;
			thiscos = cos(ang);
			thissin = sin(ang);
			imag += slope*(thissin-lsin)/(kay*kay);
			imag -= (yjbs[i]*thiscos-yjbs[i-1]*lcos)/kay;
			real += slope*(thiscos-lcos)/(kay*kay);
			real += (yjbs[i]*thissin-yjbs[i-1]*lsin)/kay;
			lsin = thissin;
			lcos = thiscos;
			}
		output(j, fundamental, real, imag, wspan);
		}

	fprintf(stderr,"Done.                                     \n");
	fclose(fp);
	return(0);
	}

