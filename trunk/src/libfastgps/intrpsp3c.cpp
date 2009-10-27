// intrpSP3c.cpp
// Methods for reading and interpolating SP3c precise ephemerides

#include "fastgps.h"
#include <sstream>
#include <iomanip>
#include <string>
#include "datetime.h"
#include "intrpsp3c.h"

#ifdef GUI
#include "../wxWidgets/fastgps_wxFrame.h"
extern fastgps_wxFrame *MainFrame;
#endif

//c::0512.08, SAH, Added new initializers and selectors.
//c::0512.09, SAH, Fixed indexing in getSP3PRNs() & getSvAccur().
//c::0711.21, SAH, Modified the code to read SP3-c files
//c::0806.30, SAH, Use Trig. Funcs. for better extrapolation

/*

 Author:  Steve Hilla, National Geodetic Survey, NOAA
 26 June 2008.

 References:
 -----------
 "A brief review of basic GPS orbit interpolation strategies",by
 Mark Schenewerk, GPS Solutions, Volume 6, Number 4, 2003. Pages 265-267.


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

// stylistic editing by Morgan Quigley, August 2008

using namespace std;
using namespace NGSdatetime;


//================== SP3cFile class =============================

//Constructors
SP3cFile::SP3cFile()
{

}

// Destructor
SP3cFile::~SP3cFile()
{
    fileStream.close();
}


//Initializers

void SP3cFile::setPathFilenameMode(char *inputFilePath)
{
  //   fileStream.open("igu14795_18.sp3.sp3",ios_base::in);
  fileStream.open(inputFilePath,std::ios::in);
  if( ! fileStream )
    fastgps_printf("Unable to open sp3 file.\n");
}


void SP3cFile::initHeaderInfo()
{
   int i,ir,j,k;

   formatVersion = ' ';
   modeFlag = ' ';
   fileType = ' ';
   SP3StartTime.SetYMDHMS(1999,1,1,0,0,0.0);
   SP3EndTime.SetYMDHMS(1999,1,1,0,0,0.0);
   numberSP3epochs = 0;
   dataUsed = "";
   coordFrame = "";
   orbitType = "";
   sourceAgency = "";
   timeSystem = "";
   gpsWeek = 9999;
   secsOfWeek = 0.0;
   SP3interval = 0.0;
   SP3mjd = 999999;
   SP3fmjd = 0.0;
   basePosVel = 0.0;
   baseClkClkrate = 0.0;
   numberSP3svs = 0;
   numberGoodPRNs = 0;
   numberGoodACCURs = 0;
   for( i = 0; i < MAXSVSEPH + 1; i++)
   {
     sp3PRNs[i] = "  0";
     svAccur[i] = 0;
   }

   lastEpochRead = 0;
   currEpoch = 0;
   numberSVparams = 0;

   for(i = 0; i < 9 + 1; i++ )
     for(j = 0; j < MAXPARAM + 1; j++ )
       for(k = 0; k < MAXSVSEPH + 1; k++ )
            inputValues[i][j][k] = 0.0;

   u_matrix = new double* [NUMTERMS];
   for(ir = 0; ir < NUMTERMS; ir++)
     u_matrix[ir] = new double[NUMTERMS];

   v_matrix = new double* [NUMTERMS];
   for(ir = 0; ir < NUMTERMS; ir++)
     v_matrix[ir] = new double[NUMTERMS];

   cvm_matrix = new double* [NUMTERMS];
   for(ir = 0; ir < NUMTERMS; ir++)
     cvm_matrix[ir] = new double[NUMTERMS];
}


bool SP3cFile::setLastEpochRead(unsigned long input)
{
   lastEpochRead = input;
   return true;
}

bool SP3cFile::setCurrEpoch(unsigned long input)
{
   currEpoch = input;
   return true;
}

bool SP3cFile::setNumberSVparams(unsigned short input)
{
   numberSVparams = input;
   return true;
}

bool SP3cFile::setNumberGoodPRNs(unsigned short input)
{
   numberGoodPRNs = input;
   return true;
}

bool SP3cFile::setNumberGoodACCURs(unsigned short input)
{
   numberGoodACCURs = input;
   return true;
}


//Selectors
int SP3cFile::readHeader()
{
  string      recType;
  string      inputRec;
  string      temp;
  char        nextFirstChar;
  long        year, month, day, hour, minute;
  double      sec;
  int         indexPRN = 0;
  int         indexACC = 0;
  int         i, prn, accur;
  bool        readLine13 = false;
  bool        readLine15 = false;

  prn = 0;
  initHeaderInfo();

  // position the stream at the beginning of the file
  fileStream.seekg(0);

  while( getline(fileStream, inputRec, '\n') )
  {
    recType = inputRec.substr( 0, 2 );
    if( recType[0] == '#' && recType[1] != '#' )
    {
      formatVersion = recType[1];
      if( formatVersion != 'a' && formatVersion != 'A' &&
          formatVersion != 'b' && formatVersion != 'B' &&
          formatVersion != 'c' && formatVersion != 'C' )
      {
        fastgps_printf("This program reads only SP3 files with format-mode: aP, bP, or cP ");
        return -1;
      }
      modeFlag = inputRec[2];
      if( modeFlag != 'P' && modeFlag != 'p' )
      {
        fastgps_printf("This program reads only SP3 files with format-mode: aP, bP, or cP ");
        return -1;
      }
      temp = inputRec.substr( 3, 4 );
      year = atol( temp.c_str() );
      temp = inputRec.substr( 8, 2 );
      month = atol( temp.c_str() );
      temp = inputRec.substr(11, 2 );
      day = atol( temp.c_str() );
      temp = inputRec.substr(14, 2 );
      hour = atol( temp.c_str() );

      temp = inputRec.substr(17, 2 );
      minute = atol( temp.c_str() );
      temp = inputRec.substr(20, 11 );
      sec = atof( temp.c_str() );
      SP3StartTime.SetYMDHMS(year, month, day, hour, minute, sec);
      SP3EndTime.SetYMDHMS(year, month, day, hour, minute, sec); // this will get added to later

      temp = inputRec.substr(32, 7 );
      numberSP3epochs = (unsigned long) atol( temp.c_str() );
      dataUsed = inputRec.substr(40, 5 );
      coordFrame = inputRec.substr(46, 5 );
      orbitType = inputRec.substr(52, 3 );
      sourceAgency = inputRec.substr(56, 4 );
    }
    else if( recType[0] == '#' && recType[1] == '#' )
    {
      temp = inputRec.substr(3, 4 );
      gpsWeek = (unsigned long) atol( temp.c_str() );
      temp = inputRec.substr(8, 15 );
      secsOfWeek = atof( temp.c_str() );
      temp = inputRec.substr(24, 14 );
      SP3interval = atof( temp.c_str() );
      temp = inputRec.substr(39, 5 );
      SP3mjd = atol( temp.c_str() );
      temp = inputRec.substr(45, 15 );
      SP3fmjd = atof( temp.c_str() );
    }
    else if( recType[0] == '+' && recType[1] == ' ' )
    {
      if( numberSP3svs == 0 )
      {
        temp = inputRec.substr(4, 2 );
        numberSP3svs = atoi( temp.c_str() );
        for( i = 9; i < 60; i=i+3 )
        {
          temp = inputRec.substr( i, 3 );
          prn = atoi( temp.substr(1,2).c_str() );
          if( prn > 0 )
          {
            sp3PRNs[indexPRN + 1] = temp;
            indexPRN++;
          }
        }
      }
      else
      {
        for( i = 9; i < 60; i=i+3 )
        {
          temp = inputRec.substr( i, 3 );
          prn = atoi( temp.substr(1,2).c_str() );
          if( prn > 0 )
          {
            sp3PRNs[indexPRN + 1] = temp;
            indexPRN++;
          }
        }
      }
    }
    else if( recType[0] == '+' && recType[1] == '+' )
    {
      for( i = 9; i < 60; i=i+3 )
      {
        temp = inputRec.substr( i, 3 );
        accur = atoi( temp.c_str() );
        if( atoi( sp3PRNs[indexACC + 1].substr(1,2).c_str() ) > 0 )
        {
          svAccur[indexACC + 1] = (unsigned short) accur;
          if( accur <= 0 )
            fastgps_printf("Accuracy warning in SP3 file processing");
          indexACC++;
        }
      }
    }
    else if( recType[0] == '%' && recType[1] == 'c' && !readLine13 )
    {
      fileType = inputRec[3];
      timeSystem = inputRec.substr( 9, 3 );
      readLine13 = true;
    }
    else if( recType[0] == '%' && recType[1] == 'f' && !readLine15 )
    {
      temp = inputRec.substr(3, 10 );
      basePosVel = atof( temp.c_str() );
      temp = inputRec.substr(14, 12 );
      baseClkClkrate = atof( temp.c_str() );
      readLine15 = true;
    }
    nextFirstChar = fileStream.peek();
    if( nextFirstChar == '*' ) break;  // exit while loop
  } // end of while loop to read all SP3 header records

  numberSVparams = 4;
  numberGoodPRNs = (unsigned short) indexPRN;
  numberGoodACCURs = (unsigned short) indexPRN;
  SP3EndTime = SP3EndTime + (double)(numberSP3epochs - 1)*SP3interval/86400.0;

  if( numberGoodPRNs <= 0  ||  numberSP3svs != numberGoodPRNs )
  {
    fastgps_printf("Error reading the PRNs from the header for SP3 file");
    return 1;
  }
  else
    return 0;
} // end of method SP3cFile::readHeader()

//------------------------------------------------------------------------

int SP3cFile::getSVPosVel(DateTime tuser, string PRNid, double rvec[])
{
  double    trun;
  double    x_new, x_dot, y_new, y_dot, z_new, z_dot;
  double    clk_new, clk_dot;
  int       i, j, k, l, m, jsv, ioverlap;
  long      jmax, jmin;
  bool      foundSV = false;
  string    inputRec;
  string    temp;
  double    SVsig, tnew, sec;
  long      year, month, day, hour, minute;
  int  iret;

  for( i = 0; i < (4 + numberSVparams); i++)
    rvec[i] = 0.0;

  for( i = 1; i <= numberSP3svs; i++)
  {
    if( PRNid == sp3PRNs[i] )
    {
      jsv = i;
      foundSV = true;
      break;  // exit loop over SP3svs once a match is found
    }
  }

  if( !foundSV )
  {
    // XPRINTF("ERROR. Cannot find PRN in SP3 file ");
    return( -2 );
  }

  SVsig = pow(2.0, (int)svAccur[jsv])/1000.0; // units: meters

  trun = (tuser - SP3StartTime)*86400.0;  // diff between DateTimes = days
  trun = trun/SP3interval + 1.0;
  jmin = (long)trun - 4;
  jmax = (long)trun + 4;

  if( jmin < 1 )
  {
    jmin = 1;
    jmax = jmin + 8;
    if( jmax > (long) numberSP3epochs)
      return( -3 );  // not enough data at beginning of SP3 file
  }

  if( jmax > (long) numberSP3epochs )
  {
    jmax = numberSP3epochs;
    jmin = jmax - 8;
    if( jmin < 1 )
      return( -4 );  // not enough data at end of SP3 file
  }

  // move in file
  if( jmax != (long)lastEpochRead   &&
      fabs( (double)lastEpochRead - (trun + 4.0) ) > (1.0/SP3interval) )
  {
    // move in file backward
    if( jmax < (long) lastEpochRead )
    {
      iret = readHeader();
      if ( iret != 0 ) 
        return iret;
    }

    // Re-use some of the data already read in, shift it to top of
    // matrix inputVlues[9+1][4+1][55+1], do the same with timeTags[] array.
    ioverlap = (int)lastEpochRead - (int)jmin + 1;
    if( (lastEpochRead > 0)  &&  (ioverlap >= 1 ) )
    {
      for( k = 1; k <= numberSP3svs; k++ )
        for( j = 1; j <= numberSVparams; j++ )
          for( i = 1; i <= (int)lastEpochRead - jmin + 1; i++ )
          {
            m = i + 9 - (lastEpochRead - jmin + 1);
            inputValues[i][j][k] = inputValues[m][j][k];
          }
      for( i = 1; i <= (int)lastEpochRead - jmin + 1; i++ )
      {
        m = i + 9 - (lastEpochRead - jmin + 1);
        timeTags[i] = timeTags[m];
      }
    } // do this section only if lastEpochRead != 0

    // skip over lines to get to the correct epoch in the SP3 file
    for( i = lastEpochRead + 1; i <= jmin - 1; i++ )
    {
      for( j = 1; j <= numberSP3svs + 1; j++ )
      {
        if( !getline(fileStream, inputRec, '\n')  )
        {
          //XPRINTF("Error skipping lines in the SP3 file.");
          return( -1 );
        }
      }
    }

    // read the data into the input values
    if( (long)(lastEpochRead + 1) > jmin)
      m = lastEpochRead + 1;
    else
      m = jmin;

    for( i = m; i <= jmax; i++ )
    {
      l = i - jmin + 1;
      if( !getline(fileStream, inputRec, '\n')  )
      {
        //XPRINTF("Error skipping time tag line in the SP3 file.");
        return( -1 );
      }
      temp = inputRec.substr( 3, 4);
      year = atol( temp.c_str() );
      temp = inputRec.substr( 8, 2);
      month = atol( temp.c_str() );
      temp = inputRec.substr(11, 2);
      day = atol( temp.c_str() );
      temp = inputRec.substr(14, 2);
      hour = atol( temp.c_str() );
      temp = inputRec.substr(17, 2);
      minute = atol( temp.c_str() );
      temp = inputRec.substr(20, 11);
      sec = atof( temp.c_str() );

      timeTags[l].SetYMDHMS(year, month, day, hour, minute, sec);

      for( k = 1; k <= numberSP3svs; k++ )
      {
        if( !getline(fileStream, inputRec, '\n')  )
        {
          //  XPRINTF("Error reading input values in the SP3 file.");
          return( -1 );
        }
        for( j = 1; j <= numberSVparams; j++ )
        {
          temp = inputRec.substr( (4 + ((j - 1)*14)), 14);
          inputValues[l][j][k] = atof( temp.c_str() );
        }
      }
    }
    lastEpochRead = jmax;

    // determine the coeffs (x_coeff, y_coeff, etc.)

    //Load the nine data values for satellite jsv
    // note that deltat is "days since deltat[0]"
    for( i = 0; i < 9; i++ )
    {
      deltat[i] = timeTags[i+1] - timeTags[1];
      xdata[i] = inputValues[i+1][1][jsv];
      ydata[i] = inputValues[i+1][2][jsv];
      zdata[i] = inputValues[i+1][3][jsv];
      clkdata[i] = inputValues[i+1][4][jsv];
      sig[i] = SVsig;
    }

    x_chisq = 0.0;
    svdfit(deltat, xdata, sig, NUMTERMS, x_coef, NUMTERMS,
           u_matrix, v_matrix, w_array, NUMTERMS, NUMTERMS,
           &x_chisq, 'P');

    y_chisq = 0.0;
    svdfit(deltat, ydata, sig, NUMTERMS, y_coef, NUMTERMS,
           u_matrix, v_matrix, w_array, NUMTERMS, NUMTERMS,
           &y_chisq, 'P');

    z_chisq = 0.0;
    svdfit(deltat, zdata, sig, NUMTERMS, z_coef, NUMTERMS,
           u_matrix, v_matrix, w_array, NUMTERMS, NUMTERMS,
           &z_chisq, 'P');

    clk_chisq = 0.0;
    svdfit(deltat, clkdata, sig, 9, clk_coef, 2,
           u_matrix, v_matrix, w_array, 9, 2,
           &clk_chisq, 'T');

  } // end of if test for jmax != lastEpochRead && ...

  // now do the actual interpolation
  tnew = tuser - timeTags[1];
  TrigExt(tnew,afunc,vfunc,NUMTERMS);

  x_new = x_dot = 0.0;
  for(j=0; j < NUMTERMS; ++j)
  {
    x_new = x_new + x_coef[j]*afunc[j];  // for position
    x_dot = x_dot + x_coef[j]*vfunc[j];  // for velocity
  }

  y_new = y_dot = 0.0;
  for(j=0; j < NUMTERMS; ++j)
  {
    y_new = y_new + y_coef[j]*afunc[j];  // for position
    y_dot = y_dot + y_coef[j]*vfunc[j];  // for velocity
  }

  z_new = z_dot = 0.0;
  for(j=0; j < NUMTERMS; ++j)
  {
    z_new = z_new + z_coef[j]*afunc[j];  // for position
    z_dot = z_dot + z_coef[j]*vfunc[j];  // for velocity
  }

  LinearFunc(tnew,afunc,vfunc,2);
  clk_new = clk_dot = 0.0;
  for(j=0; j < 2; ++j)  // num clock coeffs = 2 (A0 and A1)
  {
    clk_new = clk_new + clk_coef[j]*afunc[j];  // for SV clock
    clk_dot = clk_dot + clk_coef[j]*vfunc[j];  // for clock rate
  }

  rvec[0] = x_new;   // SV X,Y,Z in kilometers
  rvec[1] = y_new;
  rvec[2] = z_new;
  rvec[3] = clk_new;    // SV clock in microseconds
  rvec[4] = x_dot/86400.0;
  rvec[5] = y_dot/86400.0;  // velocity components in km/sec
  rvec[6] = z_dot/86400.0;
  rvec[7] = clk_dot/86400.0;    // SV clock rate in usec/sec

  return 0;
} // end of method SP3cFile::getSVPosVel()



char SP3cFile::getFormatVersion()    { return formatVersion; }

char SP3cFile::getFileType()    { return fileType; }

char SP3cFile::getModeFlag()    { return modeFlag; }

DateTime SP3cFile::getSP3StartTime()    { return SP3StartTime; }

DateTime SP3cFile::getSP3EndTime()    { return SP3EndTime; }

unsigned long SP3cFile::getNumberSP3epochs()    { return numberSP3epochs; }

string SP3cFile::getDataUsed()    { return dataUsed; }

string SP3cFile::getCoordFrame()    { return coordFrame; }

string SP3cFile::getOrbitType()    { return orbitType; }

string SP3cFile::getSourceAgency()    { return sourceAgency; }

string SP3cFile::getTimeSystem()    { return timeSystem; }

unsigned long SP3cFile::getGpsWeek()    { return gpsWeek; }

double SP3cFile::getSecsOfWeek()    { return secsOfWeek; }

double SP3cFile::getSP3interval()    { return SP3interval; }

long SP3cFile::getSP3mjd()    { return SP3mjd; }

double SP3cFile::getSP3fmjd()    { return SP3fmjd; }

double SP3cFile::getBasePosVel()    { return basePosVel; }

double SP3cFile::getBaseClkClkrate()    { return baseClkClkrate; }

unsigned short SP3cFile::getNumberSP3svs()    { return numberSP3svs; }

unsigned short SP3cFile::getSP3PRNs(string PRNids[])
{
  // Note: inside this class sp3PRNs[] is indexed starting at 1, not zero.
  // The user's array PRNids[], however, will be indexed starting at zero.
  int i;

  for( i = 1; i <=  numberSP3svs; i++ )
    PRNids[i-1] = sp3PRNs[i];

  return( numberSP3svs );

}

unsigned short SP3cFile::getSvAccur(unsigned short accurs[])
{
  // Note: inside this class svAccur[] is indexed starting at 1, not zero.
  // The user's array accurs[], however, will be indexed starting at zero.
  int i;

  for( i = 1; i <= numberSP3svs; i++ )
    accurs[i-1] = svAccur[i];

  return( numberSP3svs );

}

unsigned long SP3cFile::getLastEpochRead()    { return lastEpochRead; }

unsigned long SP3cFile::getCurrEpoch()    { return currEpoch; }

unsigned short SP3cFile::getNumberSVparams()    { return numberSVparams; }

unsigned short SP3cFile::getNumberGoodPRNs()    { return numberGoodPRNs; }

unsigned short SP3cFile::getNumberGoodACCURs()    { return numberGoodACCURs; }




void SP3cFile::TrigExt(double x, double *afunc, double *vfunc,
                       unsigned int nfunc)
{
/*
 *  Trigometric expansion basis functions.
 */

    double pi = 3.14159265358979;
    double two_pi = 2.0 * pi;

    /*
     *  Sidereal day taken from the Astronomical Almanac 1995 (USNO)
     */
    double sidereal_day = 0.99726956634;

    double period = sidereal_day;
    double P0 = two_pi / period;

    /*
     *  Reference to a specific epoch to help minimize roundoff problems.
     */
    double t0 = x;

    double P;
    int factor = 1;
    unsigned int i = 0;

    afunc[i] = 1.0; ++i;
    while( i < nfunc ) {
        P = (double)factor * P0;
        afunc[i] = sin( t0 * P ); ++i;
        afunc[i] = cos( t0 * P ); ++i;
        ++factor;
    }

    factor = 1;
    i = 0;
    vfunc[i] = 0.0; ++i;
    while( i < nfunc ) {
        P = (double)factor * P0;
        vfunc[i] = cos( t0 * P )*P; ++i;
        vfunc[i] = -1.0*sin( t0 * P )*P; ++i;
        ++factor;
    }

    return;
}

void SP3cFile::LinearFunc(double x, double *afunc, double *vfunc,
                       unsigned int nfunc)
{
/*
 *  model GPS SV clock drift with a straight line.
 *  nfunc = number of coeffs = 2
 */

    /*
     *  Reference to a specific epoch to help minimize roundoff problems.
     */
    double t0 = x;

    afunc[0] = 1.0;
    afunc[1] = t0;

    vfunc[0] = 0.0;
    vfunc[1] = 1.0;

    return;
}


bool SP3cFile::svdfit( double *X, double *Y, double *Sig,
    unsigned int NData, double *A, unsigned int MA, 
    double **U, double **V, double *W, unsigned int MP, unsigned int NP,
    double *ChiSq, char funcFlag )
{

 /*

 Modified and added to intrpsp3c.cpp in June 2008 by S. Hilla.

 Given a set of NData points X[], Y[] with individual standard
 deviations of Sig[], use chi-square minimization to determine the
 MA coefficients, A[], of the fitting function
 y = sum over i Ai * funcsi(x).
 Here we solve the fitting equation using singular value decomposition
 of the NData by MA matrix. The arrays U, V and W provide workspace
 on input. On output they define the singular value decomposition and
 can be used to obtaint he covariance matrix. MP and NP are the
 physical dimensions of the matrices U, V, and W as indicated below.
 It is necessary that MP be greater than or equal to NData and that
 NP be greather than or equal to MA. The program returns values for
 the MA fit parameters A[] and the chi-square, ChiSq. [Normally the user
 supplies a subroutine, funcs(), that returns the MA basis functions
 evaluated at x in the array afunc[], but in this version the function is
 always either TrigExt() for X,Y,Z or LinearFunc() for satellite clock]

 The input data values in *Y are checked. Bad SV clock values for a
 SP3 file are 999999 or 999999.999999, and missing epochs for a satellite's
 position may be zero filled.  If a 999999 is found for a SV clock, or
 if more than two 0.000000 values are found for a X-,Y-,Z-set of coordinates
 then the function returns false (meaning not successful).

 References:
 -----------
 "A brief review of basic GPS orbit interpolation strategies",by
 Mark Schenewerk, GPS Solutions, Volume 6, Number 4, 2003. Pages 265-267.

 Slightly modified versions of routines from
 Press, William H., Brian P. Flannery, Saul A Teukolsky and
 William T. Vetterling, 1986, "Numerical Recipes: The Art of
 Scientific Computing" (Fortran), Cambridge University Press.

 svdfit  on p. 518.
 svbksb  on pp. 57-58.
 svdcmp  on pp. 60-64.
 svdvar  on p. 519.
 */
  int sumZeroes, sum999999;
  unsigned i, j, k;
  double sum;
  double thresh;
  double tmp;
  double wmax;
  double wmin;
  double beta[nmax];
  double afunc[mmax];
  double vfunc[mmax];

  // Check the input data for svclk = 999999 and missing XYZs
  sum999999 = sumZeroes = 0;
  for( k = 0; k < NData; k++)
  {
    if( funcFlag == 'T' && (Y[k] > 999998.0) ) sum999999++;
    if( funcFlag == 'P' && (fabs(Y[k]) < 0.0000001) ) sumZeroes++;
  }
  if( sumZeroes > 1 )
  {
    cerr.setf(ios::fixed);
    cerr.setf(ios::showpoint);
    fastgps_printf("WARNING: There are data values equal to zero "); 
    for( k = 0; k < NData; k++)
    {
      cerr << k + 1 << " " << setw(20) << setprecision(15) << X[k] << "  "
        << setw(15) << setprecision(6) << Y[k] << endl;
    }
    return false;
  }

  if( sum999999 > 0 )
  {
    cerr.setf(ios::fixed);
    cerr.setf(ios::showpoint);
    fastgps_printf("WARNING: There are data values equal to 999999 "); 
    for( k = 0; k < NData; k++)
    {
      cerr << k + 1 << " " << setw(20) << setprecision(15) << X[k] << "  "
        << setw(15) << setprecision(6) << Y[k] << endl;
    }
    return false;
  }

  /* Accumulate coefficients of the fitting matrix. */
  for( i = 0; i < NData; ++i )
  {
    if( funcFlag == 'P' || funcFlag == 'p' )
      TrigExt( X[i], afunc, vfunc, MA );
    else if( funcFlag == 'T' || funcFlag == 't' )
      LinearFunc( X[i], afunc, vfunc, MA);
    else
    {
      //XPRINTF("Bad funcFlag in svdfit 1 ");
      return false;
    }

    tmp = 1.0 / Sig[i];
    for( j = 0; j < MA; ++j ) {
      U[i][j] = afunc[j] * tmp;
    }
    beta[i] = Y[i] * tmp;
  }

  /* Singular value decomposition. */
  svdcmp( U, NData, MA, MP, NP, W, V );

  /* Edit the singular values, given tol from the parameter statement,
     between here ... */
  wmax = 0.0;
  wmin = 1.0e99;
  for( j = 0; j < MA; ++j )
  {
    if( W[j] > wmax )
      wmax = W[j];
    if( W[j] < wmin )
      wmin = W[j];
  }

  thresh = tol * wmax;
  for( j = 0; j < MA; ++j ) 
    if( W[j] < thresh )
      W[j] = 0.0;

  /* ... and here. */
  svbksb( U, W, V, NData, MA, MP, NP, beta, A );

  /* Evaluate chi-square. */
  *ChiSq = 0.0;
  for( i = 0; i < NData; ++i )
  {
    if (funcFlag == 'P' || funcFlag == 'p')
      TrigExt( X[i], afunc, vfunc, MA );
    else if (funcFlag == 'T' || funcFlag == 't')
      LinearFunc( X[i], afunc, vfunc, MA );
    else
    {
      //XPRINTF("Bad funcFlag in svdfit 2 ");
      return false;
    }
    sum = 0.0;
    for( j = 0; j < MA; ++j )
      sum = sum + A[j] * afunc[j];
    tmp = ((Y[i] - sum) / Sig[i]);
    *ChiSq = *ChiSq + tmp*tmp;
  }
  return true;
}


void SP3cFile::svdvar( double **V, unsigned int MA, unsigned int NP,
    double *W, double **CVM, unsigned int NCVM )
{
  /*
     To evaluate the covariance matrix CVM of the fit for MA paramaters
     obtained by svdfit, call this routine with matrix V and W as returned
     from svdfit. NP, NCVM give the physical dimensions of V, W and CVM as
     indicated below.

     References:
     -----------
     "A brief review of basic GPS orbit interpolation strategies",by
     Mark Schenewerk, GPS Solutions, Volume 6, Number 4, 2003. Pages 265-267.
  */

  unsigned i, j, k;
  double sum;
  double wti[mmax];

  for (i = 0; i < MA; ++i)
  {
    wti[i] = 0.0;
    if( W[i] != 0.0 )
      wti[i] = 1.0 / (W[i] * W[i]);
  }

  for( i = 0; i < MA; ++i )
    for( j = 0; j <= i; ++j )
    {
      sum = 0.0;
      for( k = 0; k < MA; ++k )
        sum = sum + V[i][k] * V[j][k] * wti[k];
      CVM[i][j] = sum;
      CVM[j][i] = sum;
    }
}


void SP3cFile::svbksb( double **U, double *W, double **V, unsigned int M,
    unsigned int N, unsigned int MP, unsigned int NP,
    double *B, double *X )
{
  /*
     Solves A * X = B for a vector X where A is specified by the arrays
     U, W and V as returned by svdcmp. M and N are the logical dimensions
     of A and will be equal for a square matrices. MP and NP are the
     physical dimensions of A. B is the input right-hand side. X is the
     output solution vector. No input quantities are destroyed, so the
     routine may be called sequentially with different B's. M must be
     greater to N (see svdcmp).

     References:
     -----------
     "A brief review of basic GPS orbit interpolation strategies",by
      Mark Schenewerk, GPS Solutions, Volume 6, Number 4, 2003. Pages 265-267.
  */

  unsigned i, j;
  double S;
  double tmp[nmax];
  /* Calculate transpose U * B */
  for (j = 0; j < N; ++j)
  {
    S = 0.0;
    /* Nonzero result only if W[j] is nonzero. */
    if( W[j] != 0.0 )
    {
      for( i = 0; i < M; ++i )
        S = S + U[i][j] * B[i];
      S = S / W[j];
    }
    tmp[j] = S;
  }
  /* Multiply by V to get answer. */
  for( j = 0; j < N; ++j )
  {
    S = 0.0;
    for( i = 0; i < N; ++i )
      S = S + V[j][i] * tmp[i];
    X[j] = S;
  }
  return;
}

void SP3cFile::svdcmp( double **A, unsigned int M, unsigned int N,
    unsigned int MP, unsigned int NP, double *W, double **V )
{
  /*
     Give a matrix A, with logical dimensions M by N and physical
     dimensions MP by NP, this routine computes its singular value
     decomposition, A = U * W * transpose V. The matrix U replaces
     A on output. The diagonal matrix of singular values, W, is output
     as a vector W. The matrix V (not the transpose of V) is output as
     V. M must be greater or equal to N. If it is smaller then A should
     be filled up to square with zero rows.

     References:
     -----------
     "A brief review of basic GPS orbit interpolation strategies",by
     Mark Schenewerk, GPS Solutions, Volume 6, Number 4, 2003. Pages 265-267.
  */
  double rv1[nmax];

  /* Householder reduction to bidiagonal form. */
  int NM;
  double C, F, G = 0, H, S, X, Y, Z, Scale = 0, ANorm = 0, tmp;
  int flag;
  //unsigned i, its, j, jj, k, l;

  if (M < N)
  {
    // XPRINTF("You must augment A with extra zero rows.\n" );
    return;
  }

  for (unsigned i = 0; i < N; ++i)
  {
    unsigned l = i + 1;
    rv1[i] = Scale * G;
    G = 0.0;
    S = 0.0;
    Scale = 0.0;
    if (i < M)
    {
      for (unsigned k = i; k < M; ++k)
        Scale = Scale + fabs(A[k][i]);
      if (Scale != 0.0)
      {
        for(unsigned k = i; k < M; ++k) 
        {
          A[k][i] = A[k][i] / Scale;
          S = S + A[k][i] * A[k][i];
        }
        F = A[i][i];
        G = sqrt(S);
        if (F > 0.0)
          G = -G;
        H = F * G - S;
        A[i][i] = F - G;
        if (i != (N-1))
        {
          for (unsigned j = l; j < N; ++j)
          {
            S = 0.0;
            for (unsigned k = i; k < M; ++k)
              S = S + A[k][i] * A[k][j];
            F = S / H;
            for (unsigned k = i; k < M; ++k)
              A[k][j] = A[k][j] + F * A[k][i];
          }
        }
        for (unsigned k = i; k < M; ++k)
          A[k][i] = Scale * A[k][i];
      }
    }

    W[i] = Scale * G;
    G = 0.0;
    S = 0.0;
    Scale = 0.0;
    if (i < M && i != (N-1))
    {
      for(unsigned k = l; k < N; ++k)
        Scale = Scale + fabs(A[i][k]);
      if (Scale != 0.0)
      {
        for(unsigned k = l; k < N; ++k)
        {
          A[i][k] = A[i][k] / Scale;
          S = S + A[i][k] * A[i][k];
        }
        F = A[i][l];
        G = sqrt(S);
        if( F > 0.0 )
          G = -G;
        H = F * G - S;
        A[i][l] = F - G;
        for (unsigned k = l; k < N; ++k)
          rv1[k] = A[i][k] / H;
        if (i != (M-1))
        {
          for (unsigned j = l; j < M; ++j)
          {
            S = 0.0;
            for(unsigned k = l; k < N; ++k)
              S = S + A[j][k] * A[i][k];
            for(unsigned k = l; k < N; ++k)
              A[j][k] = A[j][k] + S * rv1[k];
          }
        }
        for (unsigned k = l; k < N; ++k)
          A[i][k] = Scale * A[i][k];
      }
    }
    tmp = fabs( W[i] ) + fabs( rv1[i] );
    if (tmp > ANorm)
      ANorm = tmp;
  }

  /* Accumulation of right-hand transformations. */
  for (int i = (int)N - 1; i >= 0; --i)
  {
    unsigned l = (unsigned)i + 1;
    if (i < (int)N - 1)
    {
      if (G != 0.0)
      {
        for (unsigned j = l; j < N; ++j)
          V[j][i] = (A[i][j] / A[i][l]) / G;
        for (unsigned j = l; j < N; ++j)
        {
          S = 0.0;
          for (unsigned k = l; k < N; ++k)
            S = S + A[i][k] * V[k][j];
          for (unsigned k = l; k < N; ++k)
            V[k][j] = V[k][j] + S * V[k][i];
        }
      }
      for (unsigned j = l; j < N; ++j)
      {
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    }
    V[i][i] = 1.0;
    G = rv1[i];
  }

  /* Accumulation of left-hand transformations. */
  for (int i = (int)N - 1; i >= 0; --i)
  {
    unsigned l = (unsigned)i + 1;
    G = W[i];
    if (i < (int)N - 1)
    {
      for (unsigned j = l; j < N; ++j)
        A[i][j] = 0.0;
    }
    if( G != 0.0 )
    {
      G = 1.0 / G;
      if (i != (int)N - 1)
      {
        for (unsigned j = l; j < N; ++j)
        {
          S = 0.0;
          for (unsigned k = l; k < M; ++k)
            S = S + A[k][i] * A[k][j];
          F = (S / A[i][i]) * G;
          for (unsigned k = i; k < M; ++k)
            A[k][j] = A[k][j] + F * A[k][i];
        }
      }
      for (unsigned j = i; j < M; ++j)
        A[j][i] = A[j][i] * G;
    }
    else
    {
      for (int j = i; j < (int)M; ++j)
        A[j][i] = 0.0;
    }
    A[i][i] = A[i][i] + 1.0;
  }

  /* Diagonalization of the bidiagonal form.
     Loop over singular values. */
  for(int k = (int)N - 1; k >= 0; --k)
  {
    /* Loop over allowed iterations. */
    for(unsigned its = 1; its <= 30; ++its)
    {
      /* Test for splitting.
         Note that rv1[0] is always zero. */
      flag = true;
      int l = k;
      for(; l >= 0; --l )
      {
        NM = l - 1;
        if ((fabs(rv1[l]) + ANorm) == ANorm)
        {
          flag = false;
          break;
        }
        else
          if ((fabs(W[NM]) + ANorm) == ANorm)
            break;
      }
      /* Cancellation of rv1[l], if l > 0; */
      if (flag)
      {
        C = 0.0;
        S = 1.0;
        for (int i = l; i <= k; ++i)
        {
          F = S * rv1[i];
          if ((fabs(F) + ANorm) != ANorm)
          {
            G = W[i];
            H = sqrt (F * F + G * G);
            W[i] = H;
            H = 1.0 / H;
            C = G * H;
            S = -F * H;
            for(unsigned j = 0; j < M; ++j )
            {
              Y = A[j][NM];
              Z = A[j][i];
              A[j][NM] = (Y * C) + (Z * S);
              A[j][i] = -(Y * S) + (Z * C);
            }
          }
        }
      }
      Z = W[k];
      /* Convergence. */
      if (l == k)
      {
        /* Singular value is made nonnegative. */
        if( Z < 0.0 )
        {
          W[k] = -Z;
          for (unsigned j = 0; j < N; ++j)
            V[j][k] = -V[j][k];
        }
        break;
      }
      if (its >= 30)
      {
        //XPRINTF("No convergence in 30 iterations.\n" );
        return;
      }

      X = W[l];
      NM = k - 1;
      Y = W[NM];
      G = rv1[NM];
      H = rv1[k];
      F = ((Y-Z)*(Y+Z) + (G-H)*(G+H)) / (2.0*H*Y);
      G = sqrt( F * F + 1.0 );
      tmp = G;
      if( F < 0.0 )
        tmp = -tmp;
      F = ((X-Z)*(X+Z) + H*((Y/(F+tmp))-H)) / X;

      /* Next QR transformation. */
      C = 1.0;
      S = 1.0;
      for (int j = l; j <= NM; ++j)
      {
        int i = j + 1;
        G = rv1[i];
        Y = W[i];
        H = S * G;
        G = C * G;
        Z = sqrt( F * F + H * H );
        rv1[j] = Z;
        C = F / Z;
        S = H / Z;
        F = (X * C) + (G * S);
        G = -(X * S) + (G * C);
        H = Y * S;
        Y = Y * C;
        for (unsigned jj = 0; jj < N; ++jj)
        {
          X = V[jj][j];
          Z = V[jj][i];
          V[jj][j] = (X * C) + (Z * S);
          V[jj][i] = -(X * S) + (Z * C);
        }
        Z = sqrt( F * F + H * H );
        W[j] = Z;

        /* Rotation can be arbitrary if Z = 0. */
        if (Z != 0.0)
        {
          Z = 1.0 / Z;
          C = F * Z;
          S = H * Z;
        }
        F = (C * G) + (S * Y);
        X = -(S * G) + (C * Y);
        for(unsigned jj = 0; jj < M; ++jj)
        {
          Y = A[jj][j];
          Z = A[jj][i];
          A[jj][j] = (Y * C) + (Z * S);
          A[jj][i] = -(Y * S) + (Z * C);
        }
      }
      rv1[l] = 0.0;
      rv1[k] = F;
      W[k] = X;
    }
  }
}

