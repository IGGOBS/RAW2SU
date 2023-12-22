/* Copyright (c) Institute of Geology and Geophysics, CAS , 202204 	*/
/* All rights reserved.		       									*/

#include "datatype.h"


/**************************** self documentation *************************/
char *sdoc[] = {
" 									",
" RAW2SU - Transfer the IGGCAS Type OBS/GOBS RAW data to SU format (V4.3)",
" 									",
" Usage:								",
" raw2su  rawfile=filename sps=* shotfile=shot.dat >stdout 		",
"	  [optinal parameters]						",
" 									",
" Required parameters:							",
" rawfile=none    input raw file name				",
" shotfile=none   input shot file name			",
" sps=*           sampling frequency			",
" 									",
" Optional parameters:							",
" t1=-5.     su data output trace start time in seconds",
" t2=30.     su data output trace end time in seconds		",
" lat=0.     instrument latitude( unit:degree )			",
" lon=0.     instrument longitude( unit:degree )			",
" offsign=1  offset is minus while shot position in the left-side/down-side of OBS station",
"            and positive in the right-side/up-side ",
"        =0  offset is absolute value",
"        =-1 positive in right/up side & negative in left/down side",
" wdep=0.0   instrument water depth(must negative, unit: m)			",
" tdcorr=0   no time correction",
"       =1   only sampling filter delay correction",
"       =2   raw file begin time correction + filter delay correction",
"       =3   all drift correction with no resampling method",
"       =4   all drift correction with resampling method",
" obsid=1    the obs deploy sequence number within line ",
" tderr=0.   time drift err from TIMEERR.LOG file (unit:s)",
" TC=-1       time control from A20**0*.LOG(e.g. A201607.LOG )file",
" rfms=-1     the micro-seconds of input rawfile in DATAFILE.LST",
" rv=0.0     reduce velcity with km/s unit,default no reduce",
" distaz=0   output the shot.no, dist, az and baz to rawfile.distaz.dat file or not (=0)",
" outfile=fhead output SU file name header, e.g. =PG04-H28, then bhx file: PG04-H28.bhx.su",
" verbose=0  =1 for more detailed information",
"	",
" 									",
" Notes:								",
" 1) This programme just transfer the OBS raw format data into su format 	",
"    according given shot file shot-by-shot with four types of time 	",
"    correction, see more details from the < OBS time instruction >	", 
" 2) The shot file has fixed <ukooa> format,please transfer the shot ",
"    information with next format:							",
"									",
"    $2 320074506.610 2001 255529.17N1180329.28E 333.3  3.33 20.  1001	",
"    id time               position		 others...		",
" 3) The tderr value should cooresponding to the rawfile, the tderr in ",
"    TIMEERR.LOG is the final error value, which means its the time error",
"    accumulated from the first file to final file. If you merged some files",
"    from the first file then the tderr=0, else if you merged the files not ",
"    from the first one, you should use the programme get_tderr get all files",
"    cooresponding time drift error.",
" 									",
NULL};

/*-----------------------------------------------------------------------
 * @ Author: Wang yuan < ywang@mail.iggcas.ac.cn >  	
 *						
 * @ IGGCAS,Beijing				
 =====
 *						
 * v4.3: Dec 18 2023
 * note: Updating the version V4.2
 *
 *------------------------------------------------------------------------*/

/***************************** end self doc *******************************/

/* main */
int main(int argc, char **argv)
{	
 	/* Define variables to be used in the command line */ 	
	float t1;			/* start time of cutted data 				*/
 	float t2;			/* end time of data after cut == cut length 		*/
	double lat;			/* obs instrument position: latitude			*/
	double lon;			/* obs instrument position: longitude			*/
	float wdep;			/* obs instrument position: water depth			*/
	float tderr;		/* time drift error from TIMEERR.LOG file		*/
	float rv;			/* reduce velocity with km/s		*/
	int sps;			/* samplint frequency of input rawfile 			*/
	int obsid;			/* the OBS deploy sequence number within line 		*/
	int tdcorr;			/* time drift correction method choice			*/
	int offsign;		/* the offset value set option:0=aboslute;1=(-,+) */
	long TC;			/* the time control value for month LOG file		*/
	int rfms;			/* the micro-seconds of input rawfile in DATAFILE.LST */
	int distaz;			/* output distaz or not ?*/
	int verbose;		/* print the information or not */
	cwp_String rawfile;		/* input file list all raw filenames 			*/
	cwp_String shotfile;		/* shot file name 					*/	
	char distazfile [100];
	cwp_String outfile;
	
	/* Define variables for filelist & shot & log file operation			*/
	int shot_num;		/* total shot number 					*/
	long long int flen;	/* file length of rawfile with bytes 			*/
	long int npts;		/* number pontints per component			*/
	int i;				/* the tmp timer				*/
	int nc;				/* number of components in rawfile		*/
	int ns;				/* per trace samples number			*/
	int ns0;			/* the tmp ns parameter				*/
	int ntr;			/* output su data total trace number		*/
	int offset;			/* distance from source point to receiver group	*/
	int shotindex;		/* shot index 0 <= shotindex < shot_num		*/
	int suFileIndex;	/* open su file or not */
	struct tm st;		/* start time of file record from rawfile name 	*/
	struct tm sh;		/* shot time 					*/
	double lon0;		/* the UTM convertion required reference longitude */	
	short zone;			/* the UTM zone */
	double stx,sty;		/* station postion coordinate -X,Y 		*/
	double shx,shy;		/* shot position coordinat - X,Y 		*/	
	double shx0,shy0;	/* start shot position coordinat - X,Y 			*/
	double shx1,shy1;	/* end shot position coordinat - X,Y 			*/
	double prox,proy;	/* project (stx,sty) to shot line => (prox,proy)*/
	long long int fpos;	/* the start read file position			*/
	long rpclk,spclk;	/* real & standard pclk from TC */
	long rfsec;
	int sdel;		/* source delay time( unit:s )			*/
	double corrstms;	/* start file time stample in micro scends */
	double shotstms;	/* shot time stample in micro scends */
	double tdms;		/* time drift value whil shot with ms unit	*/
	double total_tdms;	/* the total time drift value with ms unit	*/
	double td_slope;	/* the change slope of time drift		*/
	float t11,t22;		/* the t1 tmp 					*/
	int dtms;			/* the sampling interval with ms unit		*/
	double difftms;		/* the time difference with ms unit		*/
	double delta,az,baz;/* the source-receiver great circle, Azimuth, Back-Azimuth */
	SHOT  shot;			/* SHOT 					*/
	short rfyear;			/* the year of rawfile  */
	FILE *fpin;			/* the input raw file stream			*/
	FILE *fpshot;		/* the shot file stream				*/
	FILE *fpdistaz;
	FILE *fpx;			/* the output segy file stream -- x-components	*/
	FILE *fpy;			/* the output segy file stream -- y-components	*/
	FILE *fpz;			/* the output segy file stream -- z-components	*/
	FILE *fpd;			/* the output segy file stream -- d-components	*/
 	double *x;			/* data cut from raw file -- x component	*/
	double *y;			/* data cut from raw file -- y component	*/
	double *z;			/* data cut from raw file -- z component	*/
	double *d;			/* data cut from raw file -- d component	*/

	/* Define Variables for resample time drift correction */
	double rdt;			/* the time interval of raw data with Sec unit	*/
	double ndt;			/* the time interval of new data with Sec unit	*/
	float tN0;			/* the time difference between first read point time and shot time */
	double *rts;		/* the raw time array for raw data		*/
	double *xraw;		/* the raw time array for one trace    		*/
	double *xnew;		/* the new time array for one trace    		*/
	double *nts;		/* the new time array for add time drift	*/

	/* Shot Line direction Check */
	int shotDirect;		/* shot line main direction: NS=1;EW=2 */
	double *dshx,*dshy;	/* difference of shx,shy */
	double mdshx,mdshy;	/* mean difference of shx, shy */
	double shxo,shyo;		/* last shx, shy */
 
	/* Define variables to be used in time information */
	time_t cal_st;			/* start computing time 			*/
	struct tm *t00;			/* start computing time 			*/
	time_t cal_et;			/* end computing time 				*/
	struct tm *t99;			/* end computing time 				*/
	clock_t cpu_t;			/* the cpu total time				*/
	progress_t bar; 		/* echo progress bar 				*/
	
	/* hook up getpars */
	initargs(argc,argv);
   	requestdoc(1);	

	/* get current time as start computing time */
	time( &cal_st );
	t00 = gmtime( &cal_st);	
	
	/* get parameters */

	/* 01 -- get the four required parameters */	
	if ( !getparstring( "rawfile",&rawfile)){
		rawfile=" ";
		err( "must given the raw file name use: rawfile=filename " );	
	}
	if ( !getparstring( "shotfile",&shotfile)){
		shotfile=" ";
		warn( "must given the shot file use: shotfile=filename " );
	    	err( "if no shot file, please use programme: raw2su1 " );	
	}
	if ( !getparint("sps",&sps)){
		err("must given the sampling frequency by <sps=..>");
	}
		
	/* 02 -- get the optional parameters */
	if ( !getparint("tdcorr",&tdcorr))		tdcorr=0;
	if ( !getparint("obsid",&obsid))		obsid=1;
	if ( !getparint("offsign",&offsign))	offsign=1;
	if ( !getparfloat( "t1",  &t1  ) ) 		t1=-5.0;
	if ( !getparfloat( "t2",  &t2  ) ) 		t2=30.0;
	if ( !getpardouble( "lat", &lat ) ) 	lat=-1234.;
	if ( !getpardouble( "lon", &lon ) ) 	lon=-1234.;
	if ( !getparfloat( "wdep", &wdep) )		wdep=0.0;
	if ( !getparfloat( "tderr",&tderr ) )   tderr=0.0;
	if ( !getparfloat( "rv",&rv ) )   		rv=0.0;
	if ( !getparlong("TC",&TC ))			TC=-1;	
	if ( !getparint("rfms",&rfms))			rfms=-1;	
	if ( !getparint("verbose",&verbose))	verbose=0;	
	if ( !getparint("distaz",&distaz))		distaz=1;
	if ( !getparstring("outfile",&outfile)) outfile=rawfile;

	if ( tdcorr == 3 || tdcorr == 4 ){
	    if (tdcorr==3 && verbose ) warn( "do clock drift correction with no resampling method");
	    if (tdcorr==4 && verbose ) warn( "do clock drift correction with resampling method,need one or more minutes...");
	   	/* check for TC & OBS */
	   	if (TC == -1) err("you choose time correction method, must given TC value by TC=");
	}

	if( rv != 0. ){
		if( lat==-1234 || lon==-1234 ){ 
			warn("you choose do velocity reduction, but not given OBS lat & lon, couldn't done!!");
			rv = 0.;
		}
		
	}
	
	
/*-----------------------------
	step 0: read all file and get useful informationes 
-------------------------------*/	
	/*get input file name including information -- start recording time */
	corrstms = get_rftime( rawfile, TC, rfms, tdcorr, tderr, sps, &rfyear, &rfsec );	
	nc = get_nc( rawfile );
	npts = get_npts( rawfile, nc );
	flen = get_flen( rawfile );
	dtms = 1000./sps;
	ndt = dtms/1000.;
	rdt = ndt;

	/* get the time drift information with given parameters */ 
	//corrstms = tdcorr_begTime(tdcorr, tderr, sps, st, stms);

	if ( TC != -1 ){
		rpclk = (long)TC / 256;
		spclk = get_splck( TC );
		if (tdcorr>2){
			rdt = (spclk*ndt)/rpclk;	
		}
		if (verbose) warn("rpclk=%ld; spclk=%ld; rdt=%12.10lf, ndt=%lf",rpclk,spclk,rdt,ndt);
	} 
	
	/* get ntr & ns for output su data */
	ns = (int)( (t2-t1)*sps*1. ) + 2;
	ns0 = ns-1;
	ntr = get_ntr( npts, ns, shotfile );
	shot_num = ntr;

	/* get station postion information */
	if( lat != -1234. && lon != -1234. ){
		zone = getUTMZone(lat,lon);
		lon0 = (double)((abs(zone)-1)*6 - 180 + 3);
		geo2xy_utm(lon,lat,lon0,&stx,&sty);
	}else{
		stx=0.;
		sty=0.;
	}

	/* malloc memory */
	x = malloc( ns * sizeof(double) );
	y = malloc( ns * sizeof(double) ); 
	z = malloc( ns * sizeof(double) );
	d = malloc( ns * sizeof(double) );

	rts = malloc( ns * sizeof(double) );
	nts = malloc( ns0 * sizeof(double) );
	xnew = malloc( ns0 * sizeof(double) );

	/*initialize the time array of raw reading data */
	for (i=0; i<ns; i++)
	{
		rts[i] = i*rdt*1000;
	}

	/* open rawfile & shot file if need */
	if ( (fpin = fopen( rawfile, "r" )) == NULL )
	{
	    fprintf( stderr, "Error:cannot open rawfile:%s\n", rawfile );
    	exit(1);
	}
	
	if ( (fpshot = fopen( shotfile, "r" )) == NULL )
	{
	    fprintf( stderr, "Error:cannot open the shot file:%s\n", shotfile );
    	exit(1);
	}

	/* get the (stx,sty) project to shot line */
	/*if( stx != 0. && sty !=0. ){
		get_shotline_point(shotfile, shot_num, lon0, &shx0, &shy0, &shx1, &shy1);
		projec_piont2line(stx,sty,shx0,shy0,shx1,shy1,&prox,&proy);
	}*/

	//err("station:(%lf,%lf);shot_start:(%lf,%lf);shot_end(%lf,%lf);pro_xy(%lf,%lf)",
	//	stx,sty,shx0,shy0,shx1,shy1,prox,proy);

	/* open distaz file */
	if ( distaz )
	{
		distazfile[0]='\0';
		strcat(distazfile,rawfile);
		strcat(distazfile,".distaz.dat");
		if ( (fpdistaz = fopen( distazfile, "w" )) == NULL )
		{
	    	fprintf( stderr, "Error:cannot open distazfile:%s\n", distazfile );
    	   	exit(1);
		}
	}

	/* check shot line main direction: NS or EW */
	shotindex = 0;
	if( stx != 0. && sty !=0. ){

		/* get difference of shx,shy */
		dshx = malloc( shot_num*sizeof(double) );
		dshy = malloc( shot_num*sizeof(double) );
		while(shotindex<shot_num){

	    	shot = getshot( fpshot, rfyear );
	    	geo2xy_utm( shot.lon, shot.lat, lon0, &shx, &shy );    

	    	if ( shotindex==0 ) {
	    		dshx[shotindex] = 0.;
	    		dshy[shotindex] = 0.;
	    	} else if ( shotindex>0 ){
	    		dshx[shotindex] = fabs((shx-shxo)/1000.);
	    		dshy[shotindex] = fabs((shy-shyo)/1000.);    		
	    	}
	    	shxo = shx;
	    	shyo = shy;

	    	shotindex++;
		}

		/* rewind shot file */
		rewind( fpshot ); 
		
		/* get mean value of dshx,dshy */
		for (i = 0; i < shot_num; i++){
			mdshx += dshx[i];
			mdshy += dshy[i];
		}
		mdshx = mdshx / shot_num;
		mdshy = mdshy / shot_num;

		/* guess the main direction */
		if ( mdshx>mdshy ) shotDirect = 2; /* EW direction */
		if ( mdshx<=mdshy ) shotDirect = 1; /* NS direction */
	}

	/* main process */	
	shotindex = 0;
	suFileIndex = 0;
	fprintf( stderr,"\n" );
	/* main process loop with shot_index with two different case */
	main_process_start:	

		if (shotindex >= shot_num) goto main_process_end;

	    /* get the shot information */
	    shot = getshot( fpshot, rfyear );
	    shot.year = rfyear;
	    shotstms = get_shtm( shot,rfsec );
	   // err("%f %f",shot.lon,shot.lat);

	    /* get the distance between obs/gobs and shot */
	    if( stx != 0. && sty !=0. ){	
	    	geo2xy_utm( shot.lon, shot.lat, lon0, &shx, &shy );    
	    	offset = get_offset( stx, sty, shx, shy );
	    	get_distaz(lat, lon, shot.lat, shot.lon, &delta, &az, &baz );

	    	if( offsign == 1 ) {
	    		if( shotDirect==2 && baz>=180 ) offset *= -1; /*EW*/
	    		if( shotDirect==1 && baz>=90 && baz<270) offset *= -1; /*NS*/
	    	} else if (offsign == -1) {
	    		if( shotDirect==2 && baz<=180 ) offset *= -1; /*EW*/
	    		if( shotDirect==1 && baz>=270 && baz<=90 ) offset *= -1; /*NS*/   		
	    	} 

	    	if ( distaz )
	    	{
	    		fprintf(fpdistaz, "%d  %d  %f %f\n",shot.id,offset,az,baz );
	    	}
	    }else{
	    	offset = 0;
	    }
	    

	    /* get the start read file position	*/	    
	    fpos = get_fpos( corrstms, shotstms, tdcorr, t1, nc, dtms, rpclk, spclk, offset, rv, rdt, &tN0 );
	    //warn("%lld %lld %lf %lf ",fpos,flen,corrstms,shotstms);	     

	    /* for shot_time < file_start_time case */
	    if ( fpos < 0 )   {   
	    	if(verbose) warn("%lf %lf",corrstms,shotstms);
	    	shotindex++; 
	    	goto main_process_start;
	    }

	    /* for shot_time > file_end_time case */
	    if ( fpos > flen ) goto main_process_end;

	    /* for not enough remain raw data case */
	    if( ( flen - fpos ) < ns*nc*databytes ){
			warn( "The input rawfile not enough for the remaining shot " );
			goto main_process_end;
	    }

	    /* open su file */
	    if ( suFileIndex == 0 )
	    {
			fpx = open_sufile( outfile, "x" );
			fpy = open_sufile( outfile, "y" );
			fpz = open_sufile( outfile, "z" );
			if( nc == 4 ) fpd = open_sufile( outfile, "d" );
			suFileIndex = 1;	    	
	    }

	
	    /* cut data according given fpos and cut length */
	    read_rawfile( fpin, fpos, nc, ns, x, y, z, d );	
	    	  // err("ko");

	    /* do time drift correction with FFT interpolation method */
	    if( tdcorr == 4 )
	    {		

	    	/* initiliza the output of cutting time array */
	    	for ( i = 0; i < ns0; ++i)
	    	{
	    		nts[i] = tN0 + i*ndt*1000.;
	    	}
	    	
	    	/* do interploation & resampling by Lagrange Method */
	    	interp1d_fft(rts, x, ns, nts, xnew, ns0);	    	    	
	    	memcpy(x, xnew, ns0*sizeof(double));

	    	interp1d_fft(rts, y, ns, nts, xnew, ns0);	    	    	
	    	memcpy(y, xnew, ns0*sizeof(double));

	    	/* debugging
	    	for(i=0;i<ns;i++){
	    		fprintf(stdout,"%f %f\n",rts[i],z[i]);
	    	} */
	    	interp1d_fft(rts, z, ns, nts, xnew, ns0);	    	    	
	    	memcpy(z, xnew, ns0*sizeof(double));
	    	/* debugging 
	    	for(i=0;i<ns0;i++){
	    		fprintf(stdout,"%f %f\n",nts[i],z[i]);
	    	}*/

	    	interp1d_fft(rts, d, ns, nts, xnew, ns0);	    	    	
	    	memcpy(d, xnew, ns0*sizeof(double));	    
	    }

	
	   /*for (int i = 0; i < ns; ++i)
	    {
	       	fprintf(stdout, "%f\n", z[i] );
	    }*/
	
	    /* output su format data */	   
	    write_sudata( shot.id, obsid, shot.sn, sh, sps, ns0, wdep, t1, sdel, offset, stx, sty, shx, shy, x, fpx );
	    write_sudata( shot.id, obsid, shot.sn, sh, sps, ns0, wdep, t1, sdel, offset, stx, sty, shx, shy, y, fpy );
	    write_sudata( shot.id, obsid, shot.sn, sh, sps, ns0, wdep, t1, sdel, offset, stx, sty, shx, shy, z, fpz );
	    write_sudata( shot.id, obsid, shot.sn, sh, sps, ns0, wdep, t1, sdel, offset, stx, sty, shx, shy, d, fpd );   
	
	    /* check terminate or goto next loop? */	
	    progress_bar ( &bar, "raw2su: ", 1, shotindex+1, shot_num, 3 );	
	    shotindex++;
	    if( shotindex < shot_num ) {
	    	goto main_process_start;
	    } else {
	        goto main_process_end;
	    }
	
	/* mian process end */
	main_process_end:
	
	/* the end */
	if ( suFileIndex == 1 )
	{		
		close_file( fpx, " bhx.su\0" );
		close_file( fpy, " bhy.su\0" );
		close_file( fpz, " bhz.su\0" );
		if(nc == 4 ) close_file( fpd, " hyd.su\0" );
	}

	if (distaz)
	{
		fclose( fpdistaz );
	}
	free( x );
	free( y );
	free( z );
	if( nc == 4 ) free( d );

	free(rts);
	free(nts);
	free(xnew);

	/* get end computing time */
	time ( &cal_et );
	t99 = localtime( &cal_et );
	cpu_t = clock( );

	/* echo computing time infomation */	
	fprintf( stderr,"\n" );
	if( cpu_t/CLOCKS_PER_SEC > 15  && verbose ){	
	    fprintf( stderr,"\n" );	
	    fprintf( stderr,"\n%s: start at: %s", argv[0], asctime(t00));
	    fprintf( stderr,"%s:   end at: %s", argv[0], asctime(t99));
	    fprintf( stderr,"%s: cpu time: %lu seconds\n", argv[0], cpu_t/CLOCKS_PER_SEC ); 
	}
   	return(CWP_Exit());
}


/*--------------------------------------------------------------------------------------
	interpolate 1D seismic trace with FFT&IFFT method
----------------------------------------------------------------------------------------
INPUT:
	t0:	raw time array
	x: raw seismic trace array
	ns0: the input t0[] & x[] points;
	t1:	output time array
	ns1: the output t1[] & y[] points;
	t1: the start cut time
OUTPUT:	
	y: interpolated seismic trace array
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; April 2022
----------------------------------------------------------------------------------------*/
void interp1d_fft(double *t0,double *x, int ns0, double *t1, double *y, int ns1)
{
	int i,j;
	int *it0,*it1,*it;
	int ns;
	int ndtms,dtms,sps;
	int maxtms,nfft,nfft0,nf,nfcut;
	float ndt,fre_nyuist,ndf,onfft;

	register float *rt;		/* real trace */
	register complex *ct;	/* complex transformed trace */

	ndtms=1;

	it0 = malloc(ns0*sizeof(int));
	it1 = malloc(ns1*sizeof(int));

	/* transfer t0 & t1 to int type */
	for(i=0;i<ns0;i++){
		it0[i]=floor(t0[i]+0.5);
	}
	for(i=0;i<ns1;i++){
		it1[i]=floor(t1[i]+0.5);
	}

	dtms = it0[1]-it0[0];//warn("%d %d %d",it0[1],it0[0],dtms);
	sps = 1000/dtms;

	/* get new dens array */
	if( it1[ns1-1] >= it0[ns0-1] ){
		maxtms = it1[ns1-1];
		//warn("the output Tmax > input Tmax");
	}else{
		maxtms = it0[ns0-1];
	}
	ns = maxtms/ndtms+1;	
	it = malloc(ns*sizeof(int));
	for (i = 0; i < ns; ++i){
		it[i]=i*ndtms;
	}

	ndt = 0.001;
	fre_nyuist = sps/2.;

	/* get fft & ifft parameters */
	nfft = npfaro(ns, LOOKFAC * ns);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)  err("Padded nt=%d--too big", nfft);

	nfft0 = npfaro(ns1, LOOKFAC * ns1);
	if (nfft0 >= SU_NFLTS || nfft >= PFA_MAX)  err("Padded nt=%d--too big", nfft);

	nf = nfft/2 + 1;
	ndf = 1.0/(nfft*ndt);
	nfcut = (int)( fre_nyuist/ndf );
	onfft = 1.0/nfft0;

	/* new array with interpolate by 0 */
	rt = ealloc1float(nfft);
	ct = ealloc1complex(nf);

	memset( rt, 0.0, nfft*sizeof(float) );
	for (i = 0; i < ns0-1; i++)	{
		rt[i*dtms] = (float)x[i];
	}
	rt[ns-1] = x[ns0-1];
	//fprintf(stdout, "%f %f\n", (ns0-1)*dt,rt[new_ns-1]);

	/* do fft */
	pfarc( 1, nfft, rt, ct);

	/* set ct[i]=0 where f>Nyquist frequency */
	for( i=nfcut-1; i<nf; i++)
	{
		ct[i].r = 0.;
		ct[i].i = 0.;
	}

	/* do ifft */
	pfacr(-1, nfft, ct, rt);	

	/* debugging */
	/* for (int i = 0; i < nfft; ++i)
	{		
		fprintf(stdout, "%f\n", rt[i]*onfft );
	} */	

	/* load back and scale for inverse fft */
	memset(y,0.0,ns1*sizeof(double));
	onfft=x[0]/rt[0];
	j=0;
	i=0;	
	while (i<ns-1){
		if(it[i]==it1[j]){
			y[j] = rt[i] * onfft;
			j++;
			i++;
		}else{
			i++;
		}
	}

	if (j<ns1-1) err("fft interp error");
	y[ns1-1] = rt[ns-1]*onfft;

	//err ("%d %d %d %d %d",j,ns,ns0,ns1,ns);
	
	return;
}


void get_distaz( double stalat, double stalon, double evtlat, double evtlon, 
	double *ddelta, double *aaz, double *bbaz )
{
	/*
c
c Subroutine to calculate the Great Circle Arc distance
c    between two sets of geographic coordinates
c
c Given:  stalat => Latitude of first point (+N, -S) in degrees
c         stalon => Longitude of first point (+E, -W) in degrees
c         evtlat => Latitude of second point
c         evtlon => Longitude of second point
c
c Returns:  delta => Great Circle Arc distance in degrees
c           az    => Azimuth from pt. 1 to pt. 2 in degrees
c           baz   => Back Azimuth from pt. 2 to pt. 1 in degrees
c
c If you are calculating station-epicenter pairs, pt. 1 is the station
c
c Equations take from Bullen, pages 154, 155
c
c T. Owens, September 19, 1991
c           Sept. 25 -- fixed az and baz calculations
c
  P. Crotwell, Setember 27, 1994
            Converted to c to fix annoying problem of fortran giving wrong
               answers if the input doesn't contain a decimal point.
*/
   double delta, az, baz;
   double scolat, slon, ecolat, elon;
   double a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk;
   double rhs1,rhs2,sph,rad,del,daz,dbaz,pi,piby2;

   
   if ((stalat == evtlat)&&(stalon == evtlon)) {
      *ddelta = 0.;
      *aaz = 0.;
      *bbaz = 0.;
      return;
   }

   pi=3.141592654;
   piby2=pi/2.0;
   rad=2.*pi/360.0;
/*
c
c scolat and ecolat are the geocentric colatitudes
c as defined by Richter (pg. 318)
c
c Earth Flattening of 1/298.257 take from Bott (pg. 3)
c
*/
   sph=1.0/298.257;

   scolat=piby2 - atan((1.-sph)*(1.-sph)*tan(stalat*rad));
   ecolat=piby2 - atan((1.-sph)*(1.-sph)*tan(evtlat*rad));
   slon=stalon*rad;
   elon=evtlon*rad;
/*
c
c  a - e are as defined by Bullen (pg. 154, Sec 10.2)
c     These are defined for the pt. 1
c
*/
   a=sin(scolat)*cos(slon);
   b=sin(scolat)*sin(slon);
   c=cos(scolat);
   d=sin(slon);
   e=-cos(slon);
   g=-c*e;
   h=c*d;
   k=-sin(scolat);
/*
c
c  aa - ee are the same as a - e, except for pt. 2
c
*/
   aa=sin(ecolat)*cos(elon);
   bb=sin(ecolat)*sin(elon);
   cc=cos(ecolat);
   dd=sin(elon);
   ee=-cos(elon);
   gg=-cc*ee;
   hh=cc*dd;
   kk=-sin(ecolat);
/*
c
c  Bullen, Sec 10.2, eqn. 4
c
*/
   del=acos(a*aa + b*bb + c*cc);
   delta=del/rad;
/*
c
c  Bullen, Sec 10.2, eqn 7 / eqn 8
c
c    pt. 1 is unprimed, so this is technically the baz
c
c  Calculate baz this way to avoid quadrant problems
c
*/
   rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc - 2.;
   rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k) - 2.;
   dbaz=atan2(rhs1,rhs2);
   if (dbaz<0.0) {
      dbaz=dbaz+2*pi;
   } 
   baz=dbaz/rad;
/*
c
c  Bullen, Sec 10.2, eqn 7 / eqn 8
c
c    pt. 2 is unprimed, so this is technically the az
c
*/
   rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c - 2.;
   rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk) - 2.;
   daz=atan2(rhs1,rhs2);
   if(daz<0.0) {
      daz=daz+2*pi;
   }
   az=daz/rad;
/*
c
c   Make sure 0.0 is always 0.0, not 360.
c
*/
   if(fabs(baz-360.) < .00001) baz=0.0;
   if(fabs(az-360.) < .00001) az=0.0;

   *ddelta = delta;
   *aaz = az;
   *bbaz = baz;

   return;
   
}

/*--------------------------------------------------------------------------------------
	Get the standard Crystal oscillator frequency for different type of OBS/GOBS
----------------------------------------------------------------------------------------
INPUT:
	OBS:   	the OBS code used to identify the pclk_standard 
	TC:    	the time control value
	sps:	the sampling frequency
OUTPUT:
	pclk0:  standard Crystal oscillator frequency
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; July 2016
----------------------------------------------------------------------------------------*/
long get_splck(long TC)
{
	int i,j,k;
	long pclk0,pclk;

	pclk = (int)(TC/256/10000);

	if( pclk >= 1228 ){
	    pclk0 = 12288000;
	}else if ( pclk >= 921 ){
		pclk0 = 9216000;
	}else if ( pclk >= 307 ){
		pclk0 = 3072000;
	}

	return pclk0;
}

/*--------------------------------------------------------------------------------------
	Get the start cut file positon
----------------------------------------------------------------------------------------
INPUT:
	st  	: start recording time of input rawfile
	stms	: micro seconds of start recording time of input rawfile
	sh  	: shot time 
	shms	: micor seconds of shot time
	td_slope: the whole time drift change grandient
	t1  	: the offset value of start cut time
	nc  	: the numbers of components in rawfile
	sps 	: the sampling frequency of rawfile
OUTPUT: 
	fpos	: the file start recording postion of given shot time
	sdel	: sdel = shot time - actual cut time; unit:ms
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
long long int get_fpos(double corrstms, double shms, int tdcorr, 
			   float t1, int nc, int dtms, long rpclk, long spclk, int offset, float rv,
			   double rdt, float *tN0 )
{
	double shotms;		/* time difference with seconds unit		*/	
	time_t shot_t;		/* shot time with time_t type 				*/
	time_t file_t;		/* raw file start time with time_t type 	*/
	long long int fpos;	/* the file postion -- file pointer postion	*/
	double tdiffms;		/* time difference between shot and start recording time with mirco seconds 	*/
	float dt;		/* samplint time interval with seconds		*/
	double tdms;		/* the time drift value 					*/
	float rvt;			/* reduce velocity transfer to time rvt=D/rv*1000 */

	long long int N0;


	if ( offset==0 || rv == 0.){
		rvt = 0.; 	
	} else {
		rvt = abs(offset) * 1.0 / rv;
	}

	tdiffms = shms - corrstms + t1 * 1000.; 

	tdiffms += (int)(rvt+0.5); 

	N0 = (int)(tdiffms/(rdt*1000.)); //warn("%d,%lf,%f,%f",N0,tdiffms,rdt,tdiffms/rdt);

	fpos = N0*nc; 
	fpos = fpos*databytes;
	//err("%lf,%lf,%f,%lld,%lld,%f",corrstms,shms,tdiffms,fpos,LLONG_MAX,rvt);

	*tN0 = tdiffms - N0*rdt*1000.;

	/*** old version
	if( tdcorr > 2 ){
		tdms = (tdiffms * (rpclk - spclk)) / spclk; //err("%lf",tdms);
		tdiffms += tdms;
	}	

	tdiffms += (int)(rvt+0.5); 

	fpos = (long)( (long)((tdiffms/dtms)+0.5) * nc * databytes );

	*sdel = (int)(tdiffms - ((int)((tdiffms/dtms)+0.5)*dtms) + 0.5); 
	//err( "%lf %d %d %d\n",tdiffms,fpos,*sdel,((int)((tdiffms/dtms)+0.5)*dtms)  );

	*/

	return fpos;
}
/*--------------------------------------------------------------------------------------
	Get the raw file name inculding information
----------------------------------------------------------------------------------------
INPUT:
	rawfile : the input rawfile name
	TC	    : the time control from LOG file
	rfms    : the micorseconds of raw file begin time from DATAFILE.LST
	tdcorr  : time drift correction method
	tderr   : time drifte error from TIMEERR.LOG
	sps     : sampling frequency
	rfyear  : the year of rawfile
OUTPUT:
	stms: start time of raw file with stamp micro seconds
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
double get_rftime( char *rawfile, long TC, int rfms, int tdcorr, float tderr, int sps, short *rfyear, long *rfsec )
{
	unsigned int date[8],msecx[3],fnlength;
	short year,mon,day,jday,hour,min,sec;
	int msec,msec0;
	struct tm st;
	double sts;
	double fir_dms;
	int tderr_ms;
	double stsms;
	double corrstms;
	double pclk;
  	int i;

  	long long n1,n2;
  	long long n;
  	int m1,m2,k,c;
  	long init_TC;
  	int pclk0;

 	char a1[64]="";
  	char a2[64]="";
	char a[64]="";

	int year0,yday0,ysec0;
	int yday1,ysec1;

	//file2time( rawfile, TC, &year, &mon, &day, &jday, &hour, &min, &sec, &msec );

	/* get the date and time information for the file name */
  	fnlength = strlen ( rawfile );
  	for (i=0;i<8;i++) date[i] = HextoInt(rawfile[fnlength-12+i]);

	year  = (date[0] << 2 | date[1] >> 2)  & 0x0f;  /* 0x0f=00001111 */
  	mon   = (date[1] << 2 | date[2] >> 2)  & 0x0f;  /* 0x0f=00001111 */
  	day   = (date[2] << 3 | date[3] >> 1)  & 0x1f;  /* 0x1f=00011111 */
  	hour  = (date[3] << 4 | date[4])       & 0x1f;  /* 0x1f=00011111 */
  	min   = (date[5] << 2 | date[6] >> 2)  & 0x3f;  /* 0x3f=00111111 */
  	sec   = (date[6] << 4 | date[7])       & 0x3f;  /* 0x3f=00111111 */

  	year  = year + 2016;

  	for (i=9;i<12;i++) msecx[i-9] = HextoInt(rawfile[fnlength-12+i]);

 	msec = (msecx[0] << 8 | msecx[1] << 4 | msecx[2]);
 	msec0 = msec;

	if ( rfms != -1 )
  	{
    	pclk0 = (int)(TC/256);
    	n1 = (long long )rfms * pclk0;
    	m1 = (int)(logf(n1*1.0)/logf(2.0));
	    n2 = (long long )msec * 4096 * 10000;
    	m2 = (int)(logf(n2*1.0)/logf(2.0));
    
    	/* n1 dicemal to binary */
    	for (c = m1; c >= 0; c--)
    	{
      		k = n1 >> c;
 
      		if (k & 1)
        		a1[c]='1';
      		else
        		a1[c]='0';    
    	}
    
    	/* n2 dicemal to binary */
    	for (c = m2; c >= 0; c--)
    	{
      		k = n2 >> c;
 	
      		if (k & 1)
        		a2[c]='1';
      		else
        		a2[c]='0';    
    	}

    	if(m2>m1) {
      		strcpy(a,a1);
      		for(c=m1+1;c<=m2;c++){
        		a[c]=a2[c];
      		}  
    	}else{
      		strcpy(a,a2);
      		for(c=m2+1;c<=m1;c++){
        		a[c]=a1[c];
      		}
    	}

    	/* binary to decimal */
    	n=0;
    	for(c = m2; c >= 0; c--){
        	k=(int)a[c]-48;
        	if(k==1)  n=n+(long long )(powf(2.0,k*1.0*c));
    	}    

    	/* get intial TC */
    	init_TC = (int)(n/10000);     	
    	msec = (double)(init_TC/(pclk0*1.0)*1000);
    
  	}else{ // not given ms
    	if( TC != -1999 ) {
    		pclk = TC/256.;
    		//fprintf(stderr, "pclk=%lf\n",pclk );
      		msec = msec*4096/pclk*1000;
      		//fprintf(stderr, "%d\n",msec );
    	}  
  	}

  	year0 = year - 1;
  	yday0 = get_yday( year0, 11, 31);
  	ysec0 = yday0*24*3600;

  	yday1 = get_yday( year, mon-1, day); 
  	ysec1 = (yday1-1)*24*3600;

  	sts = ysec1 + hour*3600. + min*60. + sec; 
  	stsms = sts*1000.; 

  	if (msec0 > 1000 )	{
  		msec0 = (int)(msec0/10+0.5);
  	}

  	if (msec > 1000 ) {
  		msec = (int)(msec/10+0.5);
  	}
	
	fir_dms = FIR_POINTS*(1000./sps); 

  	//warn("%lf,%d,%d,%lf",sts,msec0,msec,fir_dms);

	if ( tdcorr == 0 ) { /* no correction */
		corrstms = (double) (stsms + msec0);		
	}else if( tdcorr == 1){ /* only sampling filter delay correction */		
		corrstms = (double) (stsms - fir_dms + msec0);
	}else if( tdcorr == 2 ) { /* rawfile name correction + filter delay correction */
		corrstms = (double) (stsms - fir_dms + msec);
	}else if( tdcorr == 3 || tdcorr == 4){ /*add tderr correction from TIMEERR.LOG */
 		corrstms = (double) (stsms - fir_dms + tderr_ms + msec);
	}else { 
		fprintf(stderr, "tdcorr=%d, tdcorr must choose between [0 4]\n",tdcorr );
		CWP_Exit();
	}

	*rfyear = year;
	*rfsec = ysec0;
	//fprintf(stderr, "Raw file begin time after correction: %04d-%02d-%02d %02d:%02d:%02d.%03d\n", 
	//	&st->tm_year+1900,&st->tm_mon+1,&st->tm_mday,&st->tm_hour, &st->tm_min, &st->tm_sec,*stms );

	return corrstms;
	 
}

/*--------------------------------------------------------------------------------------
	Transfer the rawfile name to time information
----------------------------------------------------------------------------------------
INPUT:
	rawfile	: the input rawfile name
	TC	: the time control from LOG file
OUTPUT:
	oyear : output year
	omonth: output month
	oday  : output day
	ojday : output jday
	ohour : output hour
	osec  : output seconds
	omsec : output minisecond
----------------------------------------------------------------------------------------
Credits: Zhiwei Li; IGGCAS(wuhan); Aug 2014
----------------------------------------------------------------------------------------*/
void file2time(char *fnamein, long TC,short	*oyear, short *omonth, 
	short *oday,short*ojday, short *ohour, short *omin, short *osec, int *omsec)
{ 

  unsigned int date[8],msecx[3],fnlength;
  int i;
  int year, month, day, jday, hour, min, sec, msec;
  double pclk;
	

  fnlength = strlen(fnamein);

  /* get the date and time information for the file name */
  for (i=0;i<8;i++) date[i] = HextoInt(fnamein[fnlength-12+i]);

  year  = (date[0] << 2 | date[1] >> 2)  & 0x0f;  /* 0x0f=00001111 */
  month = (date[1] << 2 | date[2] >> 2)  & 0x0f;  /* 0x0f=00001111 */
  day   = (date[2] << 3 | date[3] >> 1)  & 0x1f;  /* 0x1f=00011111 */
  hour  = (date[3] << 4 | date[4])       & 0x1f;  /* 0x1f=00011111 */
  min   = (date[5] << 2 | date[6] >> 2)  & 0x3f;  /* 0x3f=00111111 */
  sec   = (date[6] << 4 | date[7])       & 0x3f;  /* 0x3f=00111111 */

  year  = year + 2016;

  jday  = julian_day((long)year, (long)month, (long)day);
  jday  = jday - julian_day((long)year,1,1) + 1;

  for (i=9;i<12;i++) msecx[i-9] = HextoInt(fnamein[fnlength-12+i]);

  msec = (msecx[0] << 8 | msecx[1] << 4 | msecx[2]);

  if( TC != -1 ){
  	pclk = TC/256.;
  	msec = (int)(msec*4096./pclk*1000.+0.5);
  }  

  (*oyear)  = year;
  (*omonth) = month;
  (*oday)   = day;
  (*ojday)  = jday;
  (*ohour)  = hour;
  (*omin)   = min;
  (*osec)   = sec;
  (*omsec)  = msec;

  //fprintf(stderr, "Raw file begin time no correction   : %04d-%02d-%02d %02d:%02d:%02d.%03d\n",year,month,day,hour,min,sec,msec );

  return;
}  
long julian_day(long year, long mon, long day)
{return(day-32075+1461*(year+4800-(14-mon)/12)/4+367*(mon-2+(14-mon)/12*12)/12-3
*((year+4900-(14-mon)/12)/100)/4);}
/*--------------------------------------------------------------------------------------
	Get the number of componentes from input rawfile
----------------------------------------------------------------------------------------
INPUT:
	rawfile: the input rawfile name
OUTPUT:
	nc: number of componentes
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
int get_nc( char *rawfile )
{
	unsigned int date[8],fnlength;
	int i;
	short year, nc;

	fnlength = strlen(rawfile);

	/* get the date and time information for the file name */
	for (i=0;i<8;i++) date[i] = HextoInt(rawfile[fnlength-12+i]);

	year  = (date[0] << 2 | date[1] >> 2)  & 0x0f;  /* 0x0f=00001111 */
	
  	if(date[0]>>3 == 0){
  	   nc = 4;
  	}else{
  	   nc = 3;
  	}

	return nc;
}
/*--------------------------------------------------------------------------------------
	Get the number of points per data component of rawfile
----------------------------------------------------------------------------------------
INPUT:
	rawfile: the input rawfile name
	nc : number of components
OUTPUT:
	npts: number of points per data component 
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
long int get_npts( char *rawfile, int nc )
{
	long int npts;
	long long int flen;
	FILE *fpin;
	
	if ( (fpin = fopen( rawfile, "r" )) == NULL ){
		fprintf( stderr, "Error:cannot open rawfile:%s\n", rawfile );
    	exit(0);
	}
	
	fseek( fpin, 0L, 2 );
	flen = ftell( fpin );

	/* close file */
	if( fclose(fpin) !=0 )
	{
		fprintf( stderr, "Error: could not close the file %s\n",rawfile );
		exit(1);
	}

	npts = flen/nc/databytes;
	
	return npts;
}
/*--------------------------------------------------------------------------------------
	Get the file length with bytes of given filename
----------------------------------------------------------------------------------------
INPUT:
	filename: the input file name
OUTPUT:
	flen: the length of input file with bytes
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
long long int get_flen( char *filename )
{
	long long int flen;
	FILE *fpin;
	
	if ( (fpin = fopen( filename, "r" )) == NULL ){
		fprintf( stderr, "Error:cannot open rawfile:%s\n", filename );
    		exit(0);
	}
	
	fseek( fpin, 0L, 2 );
	flen = ftell( fpin );
 	rewind( fpin ); 

	/* close file */
	if( fclose(fpin) !=0 )
	{
		fprintf( stderr, "Error: could not close the file %s\n",filename );
		exit(1);
	}
	
	return flen;
}	
/*--------------------------------------------------------------------------------------
	Get the ntr and ns of output su data
----------------------------------------------------------------------------------------
INPUT:
	npts: number of points per data component
	ns: the samplint number per trace
	shotfile: the shot file name	
OUTPUT:
	ntr: number of traces for output su data
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
int get_ntr( unsigned long int npts, int ns, char *shotfile )
{
	int ntr;	/* number of traces for output su data	*/

	if( strcmp( shotfile," " ) == 0 ) {
	    ntr = (int)( npts/ns ) + 1;
	}else{
	    ntr = file_lines( shotfile );
	}
	  
 	return ntr;  	
}
/*--------------------------------------------------------------------------------------
	Open su file which will write cutted data with su format
----------------------------------------------------------------------------------------
INPUT:
	rawfile: the input rawfile name
	s: the component name
OUTPUT: 
	fpx: ouput file stream 
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Sep 2014
----------------------------------------------------------------------------------------*/
FILE *open_sufile( char *rawfile, char *s )
{
	char sufile[100];
	FILE *fp;

	/* open fpx */
	memset( sufile, 0, sizeof(sufile) );
	strcpy( sufile, rawfile );
	if( strcmp( s, "x" ) == 0 ) strcat( sufile, ".bhx.su" );
	if( strcmp( s, "y" ) == 0 ) strcat( sufile, ".bhy.su" ); 
	if( strcmp( s, "z" ) == 0 ) strcat( sufile, ".bhz.su" ); 
	if( strcmp( s, "d" ) == 0 ) strcat( sufile, ".hyd.su" ); 

   	if( ( fp = fopen( sufile,"wb+")) == NULL ){
		fprintf( stderr, " could not open the file %s\n", sufile );
		exit(1);
	}

	return fp;
}
	
/*--------------------------------------------------------------------------------------
	Read shot file only one line from given shot file stream
----------------------------------------------------------------------------------------
INPUT:
	sf: shot file stream
	year: shot year time
OUTPUT: 
	sp: SHOT information with struct SHOT
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Sep 2014
----------------------------------------------------------------------------------------*/
SHOT getshot( FILE *sf, short year )
{	
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	int slen;
	int i,j,k,h;
	int i1,j1,k1;
	char str1[30];
	SHOT sp;

	read = getline(&line, &len, sf); 	

	slen = strlen( line );

	i = 0; j = 0; k = 0; h = 0;

	char str2[slen];	
	while(i<slen)
	{	
	    j1 = -1;
	    ////fprintf( stderr,"line = %d\n",line[i] );   
	    if( line[i] != 32 && line[i]!=9 && line[i]!='\0' && line[i]!=10)
	    {
			k = k+1;		
		////fprintf( stderr,"k=%d; ",k); 		
	    }else if( (line[i]==32 || line[i]==9 || line[i]=='\0'|| line[i]==10) && k != 0 ){	
			for( j=0; j<k; j++ ) str2[j] = line[i-k+j];
			str2[k] = '\0';   
		//fprintf( stderr,"str2= %s\n", str2);
			k = 0;
			h = h + 1;
			j1 = 0;
	    }
	    i1 = -1; 
	   	//fprintf( stderr,"h=%d s=%s\n",h,str2); 	

	    if( h == 2 && j1==0 )	/* read date & time */
	    {
	    	if(strlen(str2)==12){
	    		for( j=0; j<2; j++) { i1 = i1+1; str1[j] = str2[i1]; }
				str1[j]='\0';
				sp.jday=atoi(str1); //	fprintf( stderr,"%s,day=%d; ",str1,sp.jday);
			}else{
	    		for( j=0; j<3; j++) { i1 = i1+1; str1[j] = str2[i1]; }
				str1[j]='\0';
				sp.jday=atoi(str1);	//fprintf( stderr,"%s,day=%d; ",str1,sp.jday);				
			}

		    for( j=0; j<2; j++) { i1 = i1+1; str1[j] = str2[i1]; }
       		str1[j]='\0';
       		sp.hour=atoi(str1);	//fprintf( stderr,"hour=%d; ",sp.hour);       

       		for( j=0; j<2; j++) { i1 = i1+1; str1[j] = str2[i1]; }
       		str1[j]='\0';
       		sp.min=atoi(str1);	//fprintf( stderr,"min=%d; ",sp.min);       

       		for( j=0; j<2; j++) { i1 = i1+1; str1[j] = str2[i1]; }
       		str1[j]='\0';
       		sp.sec=atoi(str1);	//fprintf( stderr,"sec=%d; ",sp.sec);	       

       		i1 = i1 + 1;
       		for( j=0; j<3; j++) { i1 = i1+1; str1[j] = str2[i1]; }
       		str1[j]='\0';
       		sp.msec=atoi(str1);	//fprintf( stderr,"msec=%d; ",sp.msec);       

       		sp.year = year;
       		kidate(sp.year, sp.jday, &sp.mon, &sp.day);       

       		//err("%d-%d-%d %d:%d:%d.%d",sp.year,sp.mon,sp.day,sp.hour,sp.min,sp.sec,sp.msec);
       		
       		j1 = -1;
	    }
		
	    if( h==3 && j1==0 ) 	/* read id */
	    {
		for( j=0;j<strlen(str2);j++) { i1 = i1+1; str1[j] = str2[i1]; }
		str1[j]='\0';
		sp.id=atoi(str1);	//fprintf( stderr,"id=%d\n ",sp.id);
	    }

	    if( h == 4 && j1==0 )	/* read Longitude & Latitude */
	    {
		k1 = -1;
		for( j=0; j<strlen(str2); j++ )
		{ 
		     k1 = k1 + 1; 
		     if( str2[k1]=='.' ) j1 = k1;
		     if( str2[k1]=='N' || str2[k1]=='S' ) break;
		}

		sp.ns[0] = str2[k1];
		sp.ns[1] = '\0'; 		//fprintf( stderr,"ns=%s,k1=%d\n",sp.ns,k1 );
	
		strmncpy( str1, str2, k1-5, k1-1 );
		sp.lasec = atof(str1);		//fprintf( stderr,"lamsec=%f\n",sp.lasec);

		strmncpy( str1, str2, k1-7, k1-6 );
		sp.lamin = atoi(str1);		//fprintf( stderr,"lamin=%d\n",sp.lamin);

		if( k1 == 9 ) strmncpy( str1, str2, k1-9, k1-8 );
		if( k1 == 8 ) strmncpy( str1, str2, k1-8, k1-8 );
		sp.ladeg = atoi(str1);		//fprintf( stderr,"ladeg=%d\n",sp.ladeg);

		sp.lat = dms2degree( sp.ladeg, sp.lamin, sp.lasec );
		
		i1 = k1;

		for( j=k1+1; j<strlen(str2); j++)
		{ 
		     k1 = k1+1; 
		     if( str2[k1]=='.' ) { j1=k1; break; }		
		}
		
		k1 = strlen(str2)-1;
		sp.ew[0] = str2[k1];
		sp.ew[1] = '\0'; 	//fprintf( stderr,"ew=%s\n",sp.ew);		
		
		strmncpy( str1, str2, k1-5, k1-1 );
		sp.losec = atof(str1);	//fprintf( stderr,"losec=%f\n",sp.losec);

		strmncpy( str1, str2, k1-7, k1-6 );
		sp.lomin = atoi(str1);	//fprintf( stderr,"lomin=%d\n",sp.lomin);

		strmncpy( str1, str2, k1-i1-1, k1-8 );
		sp.lodeg = atoi(str1);	//fprintf( stderr,"lodeg=%d\n",sp.lodeg);

		sp.lon = dms2degree( sp.lodeg, sp.lomin, sp.losec );		
	    }

	    if( h == 5 && j1==0 )	/* read sea water depth */
	    {
		sp.depth = atof(str2);	//fprintf( stderr,"depth=%f\n",sp.depth);
	    }

	    if( h == 6 && j1==0 )	/* read sea water depth */
	    {
		sp.a = atof(str2);	//fprintf( stderr,"a=%f\n",sp.a);
	    }

	    if( h == 7 && j1==0 )	/* read sea water depth */
	    {
		sp.b = atof(str2);	//fprintf( stderr,"a=%f\n",sp.b);
	    }

	    if( h == 8 && j1==0 )	/* read sea water depth */
	    {
		sp.sn = atoi(str2);	//fprintf( stderr,"a=%d\n",sp.sn);
	    }
	    i++;
	}
	
	return sp;	
}
/*--------------------------------------------------------------------------------------
	Get the offset value with unit m from given Latitude/Longitude value
----------------------------------------------------------------------------------------
INPUT:
	stx: x -- station postion
	sty: y -- station postion
	shx: x -- shot postion
	shy: y -- shot postion
OUTPUT: 
	offset  : the offset value(unit:m) which use Mercator projection transformation
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
int get_offset( double stx, double sty, double shx, double shy )
{
	int offset;
	double distance;
	   
	distance  = sqrt( pow((shx-stx),2.0) + (pow((shy-sty),2.0)) );
	if( (distance-(int)distance)>=0.5 ){ 
	    offset = (int)distance + 1; 
	}else{ 
	    offset = (int)distance; 
  	} 

	
	return offset;	
}

/*--------------------------------------------------------------------------------------
	read raw file accorrding given start postion and read length
----------------------------------------------------------------------------------------
INPUT:
	fpin: input file stream
	fpos: the start read file postion
	nc  : the number of components
	ns  : read data points 
OUTPUT: 
	xx: x-component data
	yy: y-component data
	zz: z-component data
	dd: d-component data( hydorphone )
----------------------------------------------------------------------------------------
Credits: Qinyu YOU, Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
 void read_rawfile( FILE *fp, long fpos, int nc, int ns, double *xx, double *yy, double *zz, double *dd )
{
	unsigned int i,ns1;
	NC4    data1;
	NC3    data2;

	fseek( fp, fpos, SEEK_SET );

        if( nc == 3 ){
		for( i = 0; i < ns; i++ ){
		    fread( &data2, sizeof(data2), 1, fp );
		    xx[i] = (double)dataconvert3( data2, 1 );
		    yy[i] = (double)dataconvert3( data2, 2 );
		    zz[i] = (double)dataconvert3( data2, 3 );
		}					
	}else if( nc==4 ){
		for( i = 0; i < ns; i++ ){
	 	    fread( &data1, sizeof(data1), 1, fp );
		    xx[i] = (double)dataconvert4( data1, 1 );
		    yy[i] = (double)dataconvert4( data1, 2 );
		    zz[i] = (double)dataconvert4( data1, 3 );
		    dd[i] = (double)dataconvert4( data1, 4 );
		}
	}

	return;

}
/*--------------------------------------------------------------------------------------
	get the shot time information
----------------------------------------------------------------------------------------
INPUT:
	shot: shot infomation with SHOT struct
OUTPUT:
	sh: shot time information with tm struct
	shms: micro seconds of shot time
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
double get_shtm( SHOT shot, long rfsec )
{
	int mon,day;
	struct tm sh;
	double shs;
	double shms;

	shs =  (shot.jday-1)*24*3600. + shot.hour*3600. + shot.min*60. + shot.sec*1.;
	shms = shs * 1000. + shot.msec*1.;	
	
	return shms;
}
/*--------------------------------------------------------------------------------------
	Close su file which will has been write cutted data with su format
----------------------------------------------------------------------------------------
INPUT:
	fp: ouput file stream whicn will be closed
	filename: the filename which will be closed
OUTPUT:
	none
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Nov 2014
----------------------------------------------------------------------------------------*/
void close_file( FILE *fp, char *filename )
{

	if( fclose( fp )!=0 ){
		fprintf( stderr, "Error: could not close the file:%s\n",filename );
		exit(1);
	}

	return;
}


/*--------------------------------------------------------------------------------------
	Write data & corresponding information into su file
----------------------------------------------------------------------------------------
INPUT:
	shotindex: i-th shot
	sh	 : shot time information
	sps	 : sampling frequency
	ns	 : cut data length == total cutted points
	lat	 : instrument position -- latitude
	lon	 : instrument position -- longitude
	wdep	 : instrument position -- depth
	sdel	 : shot delay time with micro seconds
	offset	 : shot to instrument distance
	stx	 : GOBS/OBS station positon -- x coordinate
	sty	 : GOBS/OBS station postion -- y coordinate
	shx	 : shot position -- x coordinate
	shy	 : shot position -- y coordinate
	x	 : the data which will write into file stream fp
	fp	 : file stream which will write data into it
OUTPUT: 
	none
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Sep 2014
----------------------------------------------------------------------------------------*/
void write_sudata(  int shotindex, int obs_no, int shot_sn, struct tm sh, int sps, int ns,  
		float wdep, float t1, int sdel, int offset, double stx, double sty,  
		double shx, double shy, double *x, FILE *fp )
{ 
	segy tr;
	int i;

	tr.tracl = shotindex;
	tr.tracr = shot_sn;
	tr.tracf = obs_no;
	tr.fldr	 = obs_no;
	tr.dt	 = (int)( 1.0/sps*1000000. + 0.5 );
	tr.ns 	 = ns;
	tr.trid	 = 1;	
	tr.sx	 = (int)shx;
	tr.sy	 = (int)shy;
	tr.gx	 = (int)stx;
	tr.gy	 = (int)sty;
	tr.offset= offset;
	tr.gelev = wdep;
	//tr.sdel  = sdel;
	tr.year  = (short)sh.tm_year + 1900;
	tr.day   = (short)sh.tm_yday;
	tr.hour  = (short)sh.tm_hour;
	tr.minute= (short)sh.tm_min;
	tr.sec   = (short)sh.tm_sec;
	tr.delrt = (int)(t1*1000);

	for( i = 0; i < ns; i ++ ) { tr.data[i] = x[i]; }// fprintf(stdout, "%f\n", tr.data[i]); }
	//err("ok");
	fputtr( fp, &tr );

	return;
}

/*--------------------------------------------------------------------------------------
	Convert the raw file binary data to float data
----------------------------------------------------------------------------------------
INPUT:
	data1: struct of raw file data type
	index: the index-th data
OUTPUT: 
	value: converted value
----------------------------------------------------------------------------------------
Credits: Xuelin Qiu, MH Zhao and SH Xia ;
	first attempt:	2008/05/30
        modified:       2009/03/16  add SPS input from the command line
        modified:       2009/04/29  merge rawb2sac and raws2sac
        modified:       2009/06/01  output npts-1, for merge in SAC to work
	South China Sea Institute of Oceanology (SCSIO), Chinese Academy of Sciences
----------------------------------------------------------------------------------------*/
float dataconvert4( NC4 data1, int index )
{
	  float value;

	 /* convert and  write the ch1 raw data */
	if( index == 1 ){
	  threeB2long.byteval[0]= data1.lbyt1;  /* swap the bytes */
	  threeB2long.byteval[1]= data1.mbyt1;
	  threeB2long.tbytev[1]= (short)data1.hbyt1;  /*sign bit extension*/
	  value = (float)threeB2long.fbytev;
	}
	
	/* convert and  write the ch2 raw data */
	if( index == 2 ){
	  threeB2long.byteval[0]= data1.lbyt2;  /* swap the bytes */
	  threeB2long.byteval[1]= data1.mbyt2;
	  threeB2long.tbytev[1]= (short)data1.hbyt2;  /*sign bit extension*/
	  value = (float)threeB2long.fbytev;
	}
	
	/* convert and write the ch3 raw data */
	if( index == 3 ){
	  threeB2long.byteval[0]= data1.lbyt3;  /* swap the bytes */
	  threeB2long.byteval[1]= data1.mbyt3;
	  threeB2long.tbytev[1]= (short)data1.hbyt3;  /*sign bit extension*/
	  value = (float)threeB2long.fbytev;
	}
	
	/* convert and write the ch4 raw data */
	if( index == 4 ){
	  threeB2long.byteval[0]= data1.lbyt4;  /* swap the bytes */
	  threeB2long.byteval[1]= data1.mbyt4;
	  threeB2long.tbytev[1]= (short)data1.hbyt4;  /*sign bit extension*/
	  value = (float)threeB2long.fbytev;
	}

	return value;
}

/*--------------------------------------------------------------------------------------
	Convert the raw file binary data to float data
----------------------------------------------------------------------------------------
INPUT:
	data2: struct of raw file data type
	index: the index-th data
OUTPUT: 
	value: converted value
----------------------------------------------------------------------------------------
Credits: Xuelin Qiu, MH Zhao and SH Xia ;
	first attempt:	2008/05/30
        modified:       2009/03/16  add SPS input from the command line
        modified:       2009/04/29  merge rawb2sac and raws2sac
        modified:       2009/06/01  output npts-1, for merge in SAC to work
	South China Sea Institute of Oceanology (SCSIO), Chinese Academy of Sciences
----------------------------------------------------------------------------------------*/
float dataconvert3( NC3 data2, int index )
{
	  float value;

	/* convert and  write the ch1 raw data */
	if( index == 1 ){
	  threeB2long.byteval[0]= data2.lbyt1;  /* swap the bytes */
	  threeB2long.byteval[1]= data2.mbyt1;
	  threeB2long.tbytev[1]= (short)data2.hbyt1;  /*sign bit extension*/

	  value = (float)threeB2long.fbytev;
	}

	/* convert and  write the ch2 raw data */
	if( index == 2 ){
	  threeB2long.byteval[0]= data2.lbyt2;  /* swap the bytes */
	  threeB2long.byteval[1]= data2.mbyt2;
	  threeB2long.tbytev[1]= (short)data2.hbyt2;  /*sign bit extension*/

	  value = (float)threeB2long.fbytev;
	}

	/* convert and write the ch3 raw data */
	if( index == 3 ){
	  threeB2long.byteval[0]= data2.lbyt3;  /* swap the bytes */
	  threeB2long.byteval[1]= data2.mbyt3;
	  threeB2long.tbytev[1]= (short)data2.hbyt3;  /*sign bit extension*/

	  value = (float)threeB2long.fbytev;
	}

	return value;
}

/* get month & day according given year & jday */
void kidate(int year, int jday, int *month, int *day)
{
  int im;
  static long ndays[12]={31,28,31,30,31,30,31,31,30,31,30,31};

  if(year == -12345 || jday == -12345) {
    *month = -12345;
    *day   = -12345;
    return;
  }

  if( (year / 4)*4 == year) {
    ndays[1] = 29;
  } else {
    ndays[1] = 28;
  }
  
  *day = jday;
  for( *month = 1; *month <= 12; (*month)++) {
    im = *month - 1;
    if(*day <= ndays[im]) 
      return;
    *day = *day - ndays[im];
  }
  *month = -12345;
  *day   = -12345;
}
/* get the yday(0-365) */
int get_yday( int year, int mon, int day )
{

	int im;
	int jday;
  	static long ndays[12]={31,28,31,30,31,30,31,31,30,31,30,31};
	
  	if( year < 0 || mon < 0 || mon > 11 || day<0 || day >31 ) {
		jday = -1;
    		return -1;
 	}
	if( year < 1900 ) year= year+1900;
	if( (year / 4)*4 == year) {
    		ndays[1] = 29;
  	} else {
    		ndays[1] = 28;
  	}

	im = 0;
	jday = 0;
	while( im < mon ){
		jday = jday + ndays[im];
		im++;
	}

	jday = jday + day;

	return jday; 
}

/* get wday(0-7) */
int get_wday( int year, int jday )
{
	int wday = 0;
	int wday11 = 0;

	if( year < 0 || jday<0 || jday >365 ) return -1;
	if( year < 1900 ) year = year + 1900;	
	wday11 = get_Jan1_wday( year );
	wday = ( wday11 + ( jday % 7 ) ) % 7;
	
	return wday;	
} 
/* Is leap year or not? */
int isLeapYear( int y)
{
 return(y % ( y % 100 ? 4 : 400 ) ? 0 : 1);
}


/* get the wday of Jan,1st */
int get_Jan1_wday( int year)
{
 	int yearNumber = year - 1;
 	int n = 0,i=0; 	

 	if(year == 0)
 	return 1;

 	for( i=1;i<year;i++)
  		if(isLeapYear(i)==1)	n++;
 	return( 365 * yearNumber + n + 1 ) % 7;
}
/*--------------------------------------------------------------------------------------
	Transfer Hexadecimal data into decimal int type
----------------------------------------------------------------------------------------
INPUT:
	databyte: the input Hexadecimal data
OUTPUT: 
	h: the decimal data after transferred
----------------------------------------------------------------------------------------
Credits: Xuelin QIU et al;
----------------------------------------------------------------------------------------*/
unsigned int HextoInt(unsigned char datebyte)
{
   unsigned int h;   

   if( !isxdigit(datebyte) ) 
   {
	fprintf( stderr,"The input filename is not a Hex number\n");
    	exit(0);
   }

   if (datebyte < 58)  h= datebyte - 48;
   else h= toupper(datebyte) - 55;

   return(h);
}



/*--------------------------------------------------------------------------------------
	Geographic coordinate convert( dd mm' ss.xxxx'' -> x.xxxx )
----------------------------------------------------------------------------------------
INPUT:
	deg: geographic degree
	min: geographic minimute
	sec: geographic second
OUTPUT: 
	degree: geographic degree
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Sep 2014
----------------------------------------------------------------------------------------*/
float dms2degree( int deg, int min, float sec )
{
	float degree;

	degree = deg + min*1.0/60.0 + sec/3600.0;
	
	return degree;
}

/*--------------------------------------------------------------------------------------
	Copies the string src[m]~src[n] to string dest[0]~dest[n-m]
----------------------------------------------------------------------------------------
INPUT:
	src: source string 
	m  : start copy position
	n  : end copy position
OUTPUT: 
	dest: copied destination string
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Sep 2014
----------------------------------------------------------------------------------------*/	
void strmncpy( char *dest, char *src, short m, short n )  
{
	short i,j;
	
	i = 0;
	for ( j = m; j <= n; j ++ )
	{
		dest[i] = src[j];
		i++;
	}
	
	dest[i] = '\0';

	return;
}

/*--------------------------------------------------------------------------------------
	Transfer given Geographic coordinate to Cartesian XY coordinate
----------------------------------------------------------------------------------------
INPUT:
	lon: longitude degree
	lat: latitude degree
OUTPUT: 
	x: x-coordinate
	y: y-coordinate
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Oct 2017
----------------------------------------------------------------------------------------*/	
void geo2xy_utm( double lon, double lat, double lon0, double *x, double *y)
{
	short zone;
	double a, f,  xx, yy;
	double xoff, ysoff, ynoff;

    a = 6378137.;
    f = 1./298.257223563;

    /* Transverse Mercator projection */
    convLLtoTM(a, f, lat, lon, lon0, &xx, &yy);

    /* Remove false Easting and Northing */
    xoff = 500000.0;
    ynoff = 0.;
    ysoff = 10000000.0;
    xx += xoff;
    yy += (zone<0) ? ysoff : ynoff;


    /* get x & y */
    *x = xx;
    *y = yy;

    return;
}

short getUTMZone(double lat, double lon)
/************************************************************************
getUTMZone - get the UTM zone number

*************************************************************************
Input:
lat     geographical latitude in degrees
lon     geographical longitude in degrees

Output:
        returns the UTM zone number

*************************************************************************
Notes:
Does computations as doubles. The latitude is positive on the northern
hemisphere and negative on the southern hemisphere. UTM coordinates
are defined between 80S and 84N. 
The longitude is negative west of the zero-meridian (Greenwich), i.e. 
its range of values is -180.0 ... 179.99999.
*************************************************************************
Author: Nils Maercklin, 30 March 2007
*************************************************************************/
{
    short zone; /* UTM zone number */
    
    /* Make sure the longitude is between -180 and 179.999 deg */
    lon = (lon+180.)-floor((lon+180.)/360.)*360.-180.;

    /* Zone number */
    zone = (short)((lon + 180.)/6.) + 1;
    if (lat >= 56.0 && lat < 64.0 && lon >= 3.0 && lon < 12.0) zone = 32;

    /* Svalbard zones */
    if (lat >= 72.0 && lat < 84.0) {
        if      (lon >= 0.0  && lon <  9.0) zone = 31;
        else if (lon >= 9.0  && lon < 21.0) zone = 33;
        else if (lon >= 21.0 && lon < 33.0) zone = 35;
        else if (lon >= 33.0 && lon < 42.0) zone = 37;
    }

    /* Return zone number */
    return zone;
}

void convLLtoTM(double a, double f, double lat, double lon, \
    double lon0, double *x, double *y)
/************************************************************************
convLLtoTM - convert latitude and longitude to Transverse Mercator 
             Easting and Northing (e.g. UTM)
*************************************************************************
Input:
a       ellipsoid simimajor axis (in m)
f       ellipsoid flattening
lat     geographical latitude in degrees
lon     geographical longitude in degrees
lon0    central meridian (longitude) in degrees

Output:
x       Easting (Transverse Mercator grid)
y       Northing (Transverse Mercator grid)

*************************************************************************
Notes:
Does computations as doubles. The longitude is negative west of the 
zero-meridian (Greenwich), i.e. its range of values is -180.0 ... 179.99999.
The latitude is positive on the northern hemisphere and negative on the 
southern hemisphere. 
For standard UTM coordinates add always a false Easting of 500 km and a 
false Northing of 10000 km on the southern hemisphere. UTM coordinates 
are only defined between 80S and 84N. 

Reference:
J. P. Snyder (1987). Map Projections - A Working Manual. 
    U.S. Geological Survey Professional Paper 1395, 383 pages.
    U.S. Government Printing Office.

This function is adopted from a Perl routine in Geo-Coordinates-UTM-0.06
by G. Crookham (CPAN, March 2007).
*************************************************************************
Author: Nils Maercklin, 30 March 2007
*************************************************************************/
{
    double e2;   /* eccentricity squared, e2 = f(2-f) */
    double ep2;  /* eccentricity prime squared */
    double k0, cn, ct, cc, ca, cm;

    /* Convert latitude and longitudes to radians */
    lat  *= PI / 180.0;
    lon  *= PI / 180.0;
    lon0 *= PI / 180.0;

    /* Ellipsoid parameters */
    e2  = f * (2.0 - f);
    ep2 = e2 / (1.0 - e2);

    /* Some constants */
    k0 = 0.9996;
    cn = a / sqrt(1.0-e2*sin(lat)*sin(lat));
    ct = tan(lat)*tan(lat);
    cc = ep2*cos(lat)*cos(lat);
    ca = cos(lat)*(lon-lon0);

    cm = a * ((1.0 - e2/4. - 3.*e2*e2/64. - 5.*e2*e2*e2/256.)*lat \
         - (3.*e2/8. + 3.*e2*e2/32.+  45.*e2*e2*e2/1024.)*sin(2.*lat) \
         + (15.*e2*e2/256. + 45*e2*e2*e2/1024.)*sin(4.*lat) \
         - (35.*e2*e2*e2/3072.)*sin(6.*lat));

    /* Transverse Mercator Easting */
    (*x) = k0 * cn * (ca + (1.-ct+cc) * ca*ca*ca/6. \
           + (5. - 18.*ct + ct*ct + 72.*cc - 58.*ep2) * ca*ca*ca*ca*ca/120.);

    /* Transverse Mercator Northing */
    (*y) = k0 * (cm + cn*tan(lat) * (ca*ca/2. + \
           (5. - ct + 9.*cc + 4.*cc*cc) * ca*ca*ca*ca/24. \
           + (61. - 58.*ct + ct*ct + 600.*cc - 330.*ep2) * \
             ca*ca*ca*ca*ca*ca/720.));
}


/*--------------------------------------------------------------------------------------
	Get total file lines from given file name 
----------------------------------------------------------------------------------------
INPUT:
	filename: the file name
OUTPUT:
	flines:	  total file lines
----------------------------------------------------------------------------------------
Credits: Yuan WANG; IGGCAS; Sep 2014
----------------------------------------------------------------------------------------*/
int file_lines( char *filename )
{
	int flines;
	int i;
	ssize_t  read;
	size_t   len = 0;
	char *   line = NULL;	
	FILE *fp;

	/* open file */
	if( (fp = fopen( filename,"r")) == NULL )
	{
	     fprintf( stderr," file_lines: could not open the file %s\n",filename );
	     return(0);
	}

	/* get lines */
	read = 0;
	flines = 0;
	while ( (read = getline(&line, &len, fp) ) > 0 )
	{
		for( i=0;i<strlen(line);i++ )
		{
		     if( line[i]!= 9 && line[i]!=32 ) break;
		}
		if( i!= strlen(line)-1 ) flines = flines + 1; 
		memset( line, 0 ,read );
	}

	/* close file */
	if( fclose(fp) !=0 )
		fprintf( stderr, "Error: could not close the file %s\n",filename );

	return flines;	
}



/* print progress bar */
void progress_bar ( progress_t *bar, char *title, int style, int index, int total, int incre )
{
	int bar_len;
	int i;
	float part_size;

	bar_len = 100/incre;	
	if( index == 0 || index == 1 ) 
	{		
		if( style<0 || style>2 ) style = 1;
		//progress_init(&bar, "", 50, PROGRESS_NUM_STYLE);
		progress_init( bar, title, bar_len, style); 
		//progress_init(&bar, "", 50, PROGRESS_BGC_STYLE);
	}	
	part_size = total*1.0/(float)bar_len;
	i = (int)( index*1.0/part_size );
	if( index == total ) { i = bar_len; }

	progress_show( bar, i, index, total );
	
	if( index == total ) progress_destroy(bar);
	
	return;
}

void progress_init(
	progress_t *bar, char *title, int max, int style)
{
    bar->chr = '=';
    bar->title = title;
    bar->style = style;
    bar->max = max;
    bar->offset = 100 / (float)max;
    bar->pro = (char *) malloc(max+1);
    if ( style == PROGRESS_BGC_STYLE )
		memset(bar->pro, 0x00, max+1);
    else {
		memset(bar->pro, 32, max);
		memset(bar->pro+max, 0x00, 1);
    }
}

void progress_show( progress_t *bar, int val, int index, int total )
{    
    int val1 = 100;
    switch ( bar->style ) 
    {
    case PROGRESS_NUM_STYLE:
	fprintf( stderr, "\033[?25l\033[30m\033[1m%s%d%%\033[?25h\033[0m\r",
		bar->title, (int)(bar->offset * val));
	fflush( stderr );
	break;
    case PROGRESS_CHR_STYLE:
	memset(bar->pro, '=', val);
	if( index == total ) 
	{
		fprintf( stderr, "\033[?25l\033[35m\033[1m%s[%-s] %d%% %d of %d \033[?25h\033[0m\r", 
			  bar->title, bar->pro, val1, index, total );
	}
	else
	{
		fprintf( stderr, "\033[?25l\033[35m\033[1m%s[%-s] %d%% %d of %d \033[?25h\033[0m\r", 
			  bar->title, bar->pro, (int)(bar->offset * val), index, total );
	}
	fflush( stderr );
	break;
    case PROGRESS_BGC_STYLE:
	memset(bar->pro, 32, val);
	fprintf( stderr, "\033[?25l\033[35m\033[1m%s\033[41m %d%% %s\033[?25h\033[0m\r", 
		bar->title, (int)(bar->offset * val), bar->pro );
	fflush( stderr );
	break;
    }
}

//destroy the the progress bar.
void progress_destroy(progress_t *bar)
{
    free(bar->pro);
}
