// intrpSP3c.h

// Class definitions for SP3-c file objects
// and prototypes for precise orbit interpolation methods

//c::0512.08,SAH,added initializers and selectors
//c::0711.21,SAH,modify the code to read the SP3-c format
//c::0806.30,SAH,use Trig.Funcs. for better extrapolation

/* notes to users:
 This class for interpolating SP3-a, SP3-b, and SP3-c files
 uses Trigonometric functions to allow for better extrapolation
 when needed.  However, such extrapolation near the beginning or
 end of a 24-hour SP3 file can still cause satellite position
 errors from 1 decimeter to one meter when extrapolating for 15
 minutes.  The error can quickly become huge after that.  So it
 is always a good idea when processing 24-hour data sets to append
 two or three SP3 files together.  That way, the interpolation
 error will always be on the order of a few mm to a few centimeters.
 These routines model the satellite clock drift as a straight line.
 The clocks for the different GPS satellites behave quite uniquely.
 The user of this source code may want to add additional logic to
 determine which satellite clocks do not have a linear drift. Any
 user using this software should test it extensively first to
 determine its accuracy for X,Y,Z coordinates and satellite clocks.
 Please send any comments to >>steve.hilla@noaa.gov<<. The matrix
 inputValues[][][] and the array timeTags[] are indexed starting at
 one to allow for easily counting through the epochs (e.g. 1 thru 96)
 found in a typical SP3 file.  The Trigonometric functions use nine
 data points (2-hours for a typical SP3 file with 15-min interval).
 This class currently only handles SP3 files with GPS SVs only, in
 the future it will need to be updated to handle mixed constellations.
 You may see a number of warnings about comparing signed and unsigned
 integers when compiling. The getSvPosVel() function returns the SV
 X,Y,Z coordinates in kilometers, the SV clock correction in micro-
 seconds, the velocity components in kilometers/second, and the SV
 clock rate-of-change in microseconds/second.

//------ Example driver program: ------------------------

#include <iostream>
#include <iomanip>
#include <string>
#include "datetime.h"
#include "intrpsp3c.h"

using namespace std;
using namespace NGSdatetime;

int main()
{
  int ierr;
  string prnID = "G04";
  double pvVec[8];
  double DeltaMjd = 300.0/86400.0;   // 300 sec = 5 min
  DateTime currTime,stopTime;
  currTime.SetYMDHMS(2008,4,2,23,45,0.0);
  stopTime.SetYMDHMS(2008,4,3,0,0,0.0);
  SP3cFile tr15min;

  tr15min.setPathFilenameMode("igs14733.sp3c",ios_base::in);
  tr15min.readHeader();
  cout << setw(18) << setprecision(6);
  while( currTime <= stopTime )
  {
    cout << "Current Time (YMDHMS): " << currTime << endl;
    ierr = (int) tr15min.getSVPosVel(currTime,prnID,pvVec);

    cout << " SV X,Y,Z (kilometers): " << pvVec[0] << " "
    << pvVec[1]  << " " << pvVec[2] << endl;
    cout << " SV clock correction (microseconds): " << pvVec[3] << endl;

    cout << " SV Xdot,Ydot,Zdot (dm/sec): " << pvVec[4]*10000.0 << " "
    << pvVec[5]*10000.0 << " " << pvVec[6]*10000.0 << endl;
    cout << " SV clock rate-of-change (10**-4 microsec/sec): "
    << pvVec[7]*10000.0 << endl;

    currTime = currTime + DeltaMjd;
  }

  cout << "Normal Termination." << endl;
  return 0;
}
//------ end of example driver program: ------------------------

-----------------------------------------------------------------------
|                                                                     |
|                  DISCLAIMER                                         |
|                                                                     |
|   THIS PROGRAM AND SUPPORTING INFORMATION IS FURNISHED BY THE       |
| GOVERNMENT OF THE UNITED STATES OF AMERICA, AND IS ACCEPTED AND     |
| USED BY THE RECIPIENT WITH THE UNDERSTANDING THAT THE UNITED STATES |
| GOVERNMENT MAKES NO WARRANTIES, EXPRESS OR IMPLIED, CONCERNING THE  |
| ACCURACY, COMPLETENESS, RELIABILITY, OR SUITABILITY OF THIS         |
| PROGRAM, OF ITS CONSTITUENT PARTS, OR OF ANY SUPPORTING DATA.       |
|                                                                     |
|   THE GOVERNMENT OF THE UNITED STATES OF AMERICA SHALL BE UNDER NO  |
| LIABILITY WHATSOEVER RESULTING FROM ANY USE OF THIS PROGRAM.  THIS  |
| PROGRAM SHOULD NOT BE RELIED UPON AS THE SOLE BASIS FOR SOLVING A   |
| PROBLEM WHOSE INCORRECT SOLUTION COULD RESULT IN INJURY TO PERSON   |
| OR PROPERTY.                                                        |
|                                                                     |
|   THIS PROGRAM IS PROPERTY OF THE GOVERNMENT OF THE UNITED STATES   |
| OF AMERICA.  THEREFORE, THE RECIPIENT FURTHER AGREES NOT TO ASSERT  |
| PROPRIETARY RIGHTS THEREIN AND NOT TO REPRESENT THIS PROGRAM TO     |
| ANYONE AS BEING OTHER THAN A GOVERNMENT PROGRAM.                    |
|                                                                     |
-----------------------------------------------------------------------

 */

#if !defined( __SP3cFILE__ )

#define __SP3cFILE__

#include <fstream>
#include <iostream>
#include <string>
#include "datetime.h"

using namespace NGSdatetime;
using namespace std;

//========================================== constants


const unsigned short  MAXSVSEPH = 33;  // for GPS-only
const unsigned short  MAXPARAM  =  8;
const unsigned short  MAXDATA    =  1000;
const unsigned short  NUMTERMS   =  9;

# define true ((int)1)
# define false ((int)0)
# define nmax ((int)1000)
# define mmax ((int)50)
# define tol ((double)1.0e-15)


//==================== SP3cFile Class =============================

class SP3cFile
{
    public:

      //Constructors
      SP3cFile();
      SP3cFile(string pathFilename);

      //Destructor
      ~SP3cFile();

      //Initializers
      void initHeaderInfo();
	  void setPathFilenameMode(char *inputFilePath);  // STG							   

      bool setLastEpochRead(unsigned long input);
      bool setCurrEpoch(unsigned long input);
      bool setNumberSVparams(unsigned short input);
      bool setNumberGoodPRNs(unsigned short input);
      bool setNumberGoodACCURs(unsigned short input);

      //Selectors
      int readHeader();
      int getSVPosVel(DateTime tuser, string PRNid, double rvec[]);

      char getFormatVersion();
      char getModeFlag();
      char getFileType();
      DateTime getSP3StartTime();
      DateTime getSP3EndTime();
      unsigned long getNumberSP3epochs();
      string getDataUsed();
      string getCoordFrame();
      string getOrbitType();
      string getSourceAgency();
      string getTimeSystem();
      unsigned long getGpsWeek();
      double getSecsOfWeek();
      double getSP3interval();
      long getSP3mjd();
      double getSP3fmjd();
      double getBasePosVel();
      double getBaseClkClkrate();
      unsigned short getNumberSP3svs();
      unsigned short getSP3PRNs(string PRNids[MAXSVSEPH + 1]);
      unsigned short getSvAccur(unsigned short accurs[MAXSVSEPH + 1]);
      unsigned long getLastEpochRead();
      unsigned long getCurrEpoch();
      unsigned short getNumberSVparams();
      unsigned short getNumberGoodPRNs();
      unsigned short getNumberGoodACCURs();

void TrigExt(double x, double *afunc, double *vfunc, unsigned int nfunc);
void LinearFunc(double x, double *afunc, double *vfunc, unsigned int nfunc);
bool svdfit( double *X, double *Y, double *Sig, unsigned int NData,
    double *A, unsigned int MA,
    double **U, double **V, double *W, unsigned int MP, unsigned int NP,
    double *ChiSq, char funcFlag );
void svdcmp( double **A, unsigned int M, unsigned int N,
    unsigned int MP, unsigned int NP, double *W, double **V );
void svbksb( double **U, double *W, double **V, unsigned int M,
    unsigned int N, unsigned int MP, unsigned int NP,
    double *B, double *X );
void svdvar( double **V, unsigned int MA, unsigned int NP,
    double *W, double **CVM, unsigned int NCVM );

    private:
      string               pathFilename;
      fstream              fileStream;
      ios_base::openmode   fileMode;

      char                 formatVersion;
      char                 modeFlag;
      char                 fileType;    // e.g. G = GPS, M = MIXED, R = GLONASS
      DateTime             SP3StartTime;
      DateTime             SP3EndTime;
      unsigned long        numberSP3epochs;
      string               dataUsed;
      string               coordFrame;
      string               orbitType;
      string               sourceAgency;
      string               timeSystem;   // e.g. GPS,GLO,GAL,TAI,UTC
      string               sp3PRNs[MAXSVSEPH + 1];
      unsigned long        gpsWeek;
      unsigned long        count_ntrpl;
      double               secsOfWeek;
      double               SP3interval;
      long                 SP3mjd;
      double               SP3fmjd;
      double               basePosVel;        // e.g. 1.2500000
      double               baseClkClkrate;    // e.g. 1.025000000
      unsigned short       numberSP3svs;
      unsigned short       svAccur[MAXSVSEPH + 1];

      unsigned long        lastEpochRead;
      unsigned long        currEpoch;
      unsigned short       numberSVparams;
      DateTime             timeTags[9 + 1];
      double               inputValues[9 + 1][MAXPARAM + 1][MAXSVSEPH + 1];
      double               deltat[NUMTERMS];
      double               xdata[NUMTERMS];
      double               ydata[NUMTERMS];
      double               zdata[NUMTERMS];
      double               clkdata[NUMTERMS];
      double               sig[NUMTERMS];
      double               x_coef[NUMTERMS];
      double               y_coef[NUMTERMS];
      double               z_coef[NUMTERMS];
      double               clk_coef[NUMTERMS];
      double               afunc[NUMTERMS];
      double               vfunc[NUMTERMS];
      double               w_array[NUMTERMS];
      double               x_chisq, y_chisq, z_chisq, clk_chisq;
      double               **u_matrix;
      double               **v_matrix;
      double               **cvm_matrix;
      unsigned short       numberGoodPRNs;
      unsigned short       numberGoodACCURs;

};

#endif
