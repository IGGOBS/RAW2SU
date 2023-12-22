/***************************************************************************************
	Define all kinds of data type with struct or union format
****************************************************************************************
* 
****************************************************************************************
Credits: 
	Yuan WANG;
	@ywang@mail.iggcas.ac.cn
	Institute of Geology and Geophysics,Chinese Acadamy of Sciences;
	Sep 09 2014
****************************************************************************************/  
#include <time.h>
#include <stdbool.h>
#include <stdio.h>

#include "cwp.h"
#include "su.h"
#include "segy.h"
#include "progress.h"

#define NPTS 50000000
#define file_len  20			/* raw file name string length			*/
#define databytes  3			/* single record data byte number		*/
#define FIR_POINTS 18		/* fitler delay samplint points */
#define tdsign  1			/* time dirft correction sign: 
					   1: OBS time - true time > 0 <--> OBS time faster <--> pclk < pclk0
					  -1: OBS time - true time > 0 <--> OBS time slower <--> pclk > pclk0
					*/ 


#define PFA_MAX	720720	/* Largest allowed nfft	          */
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */

/* SHOT file struct */
typedef struct 
{
	int   year;
	int   mon;
	int   day;
	int   jday;
	int   hour;
	int   min;
	int   sec;
	int   msec;
	int   id;
	char  ns[2];
	float lasec;
	int   lamin;
	int   ladeg;
	char  ew[2];
	float losec;
	int   lomin;
	int   lodeg;
	float depth;
	float a;
	float b;
	int   sn;
	double lat;
	double lon;
	float reflat;
	double x;
	double y;
}SHOT;

typedef struct
{ 
	char hbyt1,mbyt1,lbyt1;		/* x value:E radial direction	*/
	char hbyt2,mbyt2,lbyt2;		/* y value:N transverse direction	*/
	char hbyt3,mbyt3,lbyt3;		/* z value:V vertical direction	*/
	char hbyt4,mbyt4,lbyt4;		/* h value:hydrophone 		*/
} NC4;

typedef struct
{
	char hbyt1,mbyt1,lbyt1; 	/* x value:E radial direction	*/
	char hbyt2,mbyt2,lbyt2;		/* y value:N transverse direction	*/
     	char hbyt3,mbyt3,lbyt3;		/* z value:V vertical direction	*/
} NC3;

union
{
	unsigned char byteval[4];
    	short tbytev[2];
     	int   fbytev; // lzw----   for 64bit system
} threeB2long;


/* Internal subroutine statment */
void close_file( FILE *fpx, char *filename );
void convLLtoTM(double a, double f, double lat, double lon, \
    double lon0, double *x, double *y);
float dataconvert3( NC3 data2, int index );
float dataconvert4( NC4 data1, int index );
float dms2degree( int deg, int min, float sec );
void file2time(char *fnamein, long TC,short	*oyear, short *omonth, 
	short *oday,short*ojday, short *ohour, short *omin, short *osec, int *omsec);
int file_lines( char *filename );
void geo2xy_utm( double lon, double lat, double lon0, double *x, double *y);
void get_distaz( double stalat, double stalon, double evtlat, double evtlon, 
	double *ddelta, double *aaz, double *bbaz );
long long int get_flen( char *filename );
long long int get_fpos( double corrstms, double shms, int tdcorr,
			   float t1, int nc, int dtms, long rpclk, long spclk, int offset,
			   float rv, double rdt, float *tN0 );
int get_Jan1_wday( int year);
int get_nc( char *rawfile );
long int get_npts( char *rawfile, int nc );
int get_ntr( unsigned long int npts, int ns, char *shotfile );
int get_offset( double stx, double sty, double shx, double shy );
double get_rftime( char *rawfile, long TC, int rfms, int tdcorr, float tderr, int sps, short *year, long *rfsec );
SHOT getshot( FILE *s, short year );
double get_shtm( SHOT shot, long rfsec);
long get_splck( long TC );
short getUTMZone(double lat, double lon);
int get_wday( int year, int jday );
int get_yday( int year, int mon, int day );
unsigned int HextoInt(unsigned char datebyte);
int ibitr( int j, int nu);
void interp1d_fft(double *t0,double *x, int ns0, double *t1, double *y, int ns1);
int isLeapYear( int y);
long julian_day(long year, long mon, long day);
time_t mktime_t( int year, int mon, int day, int hour, int min, int sec );
void kidate(int year, int jday, int *month, int *day);
struct tm mktm( int year, int mon, int day, int hour, int min, int sec );
FILE *open_sufile( char *rawfile, char *s );
void progress_bar ( progress_t *bar, char *title, int style, int index, int total, int incre );
void read_rawfile( FILE *fp, long fpos, int nc, int ns, double *xx, double *yy, double *zz, double *dd );
void resample_tdcorr( int npts, float tolms, float *xin, float *xout, float *y );
char * strmcpy( char *s, int m0, int m2 );
void strmncpy( char *dest, char *src, short m, short n );
double tdcorr_begTime(int tdcorr, float tderr, int sps, struct tm st, int stms);
double tmdiff( struct tm t1, struct tm t0 );
void write_sudata(  int shotindex, int obs_no, int shot_sn, struct tm sh, int sps, int ns,
		float wdep, float t1, int sdel, int offset, double stx, double sty, double shx, double shy, 
		double *x, FILE *fp );
void update_shottime( struct tm *sh, float t1 );
