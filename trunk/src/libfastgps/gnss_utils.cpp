/*
* Copyright (c) 2008, Scott Gleason
* All rights reserved.
*
* Originally written by Scott Gleason

* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the authors' names nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>

/****************************************************************************
* Function: void wgsllh2xyz(double *llh,double *xyz)
*
* Converts from WGS84 latitude, longitude and height into WGS X, Y and Z
* From: Simo H. Laurila, "Electronic Surveying and Navigation", John Wiley & Sons (1976).
****************************************************************************/
void wgsllh2xyz(double *llh,double *xyz)
{
  const double a = 6378137.0E0;  // semi-major axis
  const double e = 0.0818191908426E0;   // eccentricity
  const double ome2 = 0.99330562000987; // (1.0E0 - e2)
  double n,d,nph,tempd;

  d = e*sin(llh[0]);
  n = a/sqrt(1.0E0-d*d);
  nph = n + llh[2];

  tempd = nph*cos(llh[0]);
  xyz[0] = tempd*cos(llh[1]);
  xyz[1] = tempd*sin(llh[1]);
  xyz[2] = (ome2*n + llh[2])*sin(llh[0]);
}

/****************************************************************************
* Function: void wgsxyz2llh(double *xyz,double *llh)
*
* Converts from WGS84 X,Y Z to lat,lon, hgt
* From: Simo H. Laurila, "Electronic Surveying and Navigation", John Wiley & Sons (1976).
****************************************************************************/
void wgsxyz2llh(double *xyz,double *llh)
{
  const double tol = 1.0E-13;
  const double dpi2 = 1.570796326794897E0;
  const double a = 6378137.0E0;		    // semi-major axis
  const double b = 6356752.3142E0;      // semi-minor axis
  const double e = 0.0818191908426E0;   // eccentricity
  const double e2 = 0.00669437999013E0; // e^2

  double n,d,nph,tempd,latx,hgtx;

  // see if we are close to the Z axis
  tempd = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
  if (tempd <= tol)
  {
    // we are very close to the Z axis
    llh[0] = dpi2;
    if(xyz[2] < 0.0){
      llh[0] = -llh[0];
    }
    llh[1] = 0.0;
    llh[1] = fabs(xyz[2]) - b;     /* height at poles */
    return;
  }

  // for all other cases
  llh[0] = atan2(xyz[2],tempd);
  llh[1] = atan2(xyz[1],xyz[0]);

  // height needs to be iterated
  llh[2] = tempd/cos(llh[0]);
  latx = llh[0] + 1.0;
  hgtx = llh[2] + 1.0;

  while(fabs(llh[0] - latx) >= tol || fabs(llh[2] - hgtx) >= 0.01)
  {
    latx = llh[0];
    hgtx = llh[2];
    d = e*sin(latx);
    n = a/sqrt(1.0 - d*d);
    llh[2] = tempd/cos(latx) - n;
    nph = n + llh[2];
    d = 1.0 - e2*(n/nph);
    llh[0] = atan2(xyz[2],tempd*d);
  }
}

/**************************************************************************** */
//
// Converts wgs84 ecef coordinated into a local north,east,up frame
// rotated around the reference wgs84 ecef position
//
/****************************************************************************/
void wgsxyz2neu(double ecef[3],double ref_ecef[3],double NEU[3])
{
  double M[3][3];
  double ref_el,ref_az;
  double tempd;

  // convert reference point to spherical earth centered coords
  tempd = sqrt(ref_ecef[0]*ref_ecef[0] + ref_ecef[1]*ref_ecef[1]);
  ref_el = atan2(ref_ecef[2],tempd);
  ref_az = atan2(ref_ecef[1],ref_ecef[0]);

  M[0][0] = -sin(ref_el)*cos(ref_az);
  M[0][1] = -sin(ref_el)*sin(ref_az);
  M[0][2] = cos(ref_el);
  M[1][0] = -sin(ref_az);
  M[1][1] = cos(ref_az);
  M[1][2] = 0.0;
  M[2][0] = cos(ref_el)*cos(ref_az);
  M[2][1] = cos(ref_el)*sin(ref_az);
  M[2][2] = sin(ref_el);

  NEU[0] = M[0][0]*ecef[0] + M[0][1]*ecef[1] + M[0][2]*ecef[2];
  NEU[1] = M[1][0]*ecef[0] + M[1][1]*ecef[1] + M[1][2]*ecef[2];
  NEU[2] = M[2][0]*ecef[0] + M[2][1]*ecef[1] + M[2][2]*ecef[2];
}





