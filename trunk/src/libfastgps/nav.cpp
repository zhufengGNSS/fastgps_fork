/*
* Copyright (c) 2008, Morgan Quigley, Pieter Abbeel and Scott Gleason
* All rights reserved.
*
* Written by Scott Gleason
* Additional contributions by Morgan Quigley and Pieter Abbeel
*
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

#include "fastgps.h"

#ifdef GUI
#include "../wxWidgets/fastgps_wxFrame.h"
extern fastgps_wxFrame *MainFrame;
#endif


unsigned int PVT_Solution(int num_chans)
{
  unsigned int num_valid_chans = 0;
  unsigned int retval = SUCCESS;
  s_PVT_Info  PVT_Info[MAX_CHANNELS];
  double Rx_Pos[3];
  double Rx_Clock_Bias;
  unsigned int max_iterations,i,j,k;
  double omp[MAX_CHANNELS];
  double p_pred[MAX_CHANNELS];
  double G[MAX_CHANNELS][4];
  double Gtrans[4][MAX_CHANNELS];
  double GtG[4][4];
  double tempd,v0vk_mag,tau,wEtau;
  double tempv[3];
  double tempvX[MAX_CHANNELS];
  double rotm[3][3];
  double vk_new[3];
  double los[3];
  double correction[4];
  double velsol[4];
  double Vs[MAX_CHANNELS][4];
  double pdot[MAX_CHANNELS];
  double H[4][4];
  double X[4][MAX_CHANNELS];
  double temp_XYZ[3];

  correction[0] = 0.0;
  correction[1] = 0.0;
  correction[2] = 0.0;
  correction[3] = 0.0;

  // extract valid measurements at this time epoch
  for (int i=0;i<num_chans;i++)
  {
    if(c[i].nav.valid_for_pvt == NO)
      continue;

    PVT_Info[num_valid_chans].satposxyz[0] = c[i].nav.sat_pos[0];
    PVT_Info[num_valid_chans].satposxyz[1] = c[i].nav.sat_pos[1];
    PVT_Info[num_valid_chans].satposxyz[2] = c[i].nav.sat_pos[2];
    PVT_Info[num_valid_chans].satclk_bias = c[i].nav.clock_corr;
    PVT_Info[num_valid_chans].satvelxyz[0] = c[i].nav.sat_vel[0];
    PVT_Info[num_valid_chans].satvelxyz[1] = c[i].nav.sat_vel[1];
    PVT_Info[num_valid_chans].satvelxyz[2] = c[i].nav.sat_vel[2];
    PVT_Info[num_valid_chans].satclk_drift = c[i].nav.clock_drift;
    PVT_Info[num_valid_chans].pseudorange = c[i].nav.pseudorange;
    PVT_Info[num_valid_chans].doppler = c[i].nav.doppler_meas;

    // apply sat clock bias correction to psuedorange
    PVT_Info[num_valid_chans].pseudorange += 
                                   PVT_Info[num_valid_chans].satclk_bias*NAV_C;
    // calculate pseudorange rate, m/s
    PVT_Info[num_valid_chans].pseudorange_dot = 
                     ((NAV_C*PVT_Info[num_valid_chans].doppler)/(CARRIER_FREQ))
                     - PVT_Info[num_valid_chans].satclk_drift*NAV_C;
    // TODO: velocity info
    num_valid_chans++;
  }

  if(num_valid_chans < MINIMUM_PVT_SATELLITES)
    return(FAILURE);

  //####################
  // Position Solution #
  //####################

  // Initial guess for Rx position, for starters lets use the truth with an 
  // offset (so as not to make it too easy)
  Rx_Pos[0] = system_vars.initial_pos_guess[0];   // WGS84 X
  Rx_Pos[1] = system_vars.initial_pos_guess[1];   // WGS84 Y
  Rx_Pos[2] = system_vars.initial_pos_guess[2];   // WGS84 Z

  // if solution does not converge after this many iterations, 
  // something is wrong
  max_iterations = 10;

  for (i=0;i<max_iterations;i++)
  {
    // loop through valid measurements
    for (j=0;j<num_valid_chans;j++)
    {
      // The satellite positions need to be corrected for Earth Rotation 
      // during the transmission time.
      vector_subtract(Rx_Pos,PVT_Info[j].satposxyz,tempv);
      v0vk_mag = vector_norm(tempv);
      tau = v0vk_mag/NAV_C;
      wEtau = NAV_OMEGAE_DOT*tau;
      rotm[0][0] = cos(wEtau);
      rotm[0][1] = sin(wEtau);
      rotm[0][2] = 0.0;
      rotm[1][0] = -sin(wEtau);
      rotm[1][1] = cos(wEtau);
      rotm[1][2] = 0.0;
      rotm[2][0] = 0.0;
      rotm[2][1] = 0.0;
      rotm[2][2] = 1.0;
      matrix_multiply(3, 3, 1, (double *)rotm,
                      (double *)PVT_Info[j].satposxyz,
                      (double *)vk_new);  // result in vk_new
      // predicted range from satellite position and estimated Rx position
      vector_subtract(Rx_Pos,vk_new,tempv);
      p_pred[j] = vector_norm(tempv);
      // observed minus predicted range
      omp[j] = PVT_Info[j].pseudorange - p_pred[j];
      // line of sight unit vector
      vector_subtract(PVT_Info[j].satposxyz,Rx_Pos,los);
      los[0] = los[0]/p_pred[j];
      los[1] = los[1]/p_pred[j];
      los[2] = los[2]/p_pred[j];

      // build geometry matrix
      G[j][0] = -1*los[0];
      G[j][1] = -1*los[1];
      G[j][2] = -1*los[2];
      G[j][3] = 1;

    } // end of  channel loop
    // solve for position corrections
    matrix_transpose(num_valid_chans,4,(double *)G,(double *)Gtrans);
    matrix_multiply(4,num_valid_chans,4,(double *)Gtrans,(double *)G,
                    (double *)GtG);
    invert4x4(GtG,H);
    matrix_multiply(4,4,num_valid_chans,(double *)H,(double *)Gtrans,
                    (double *)X);
    matrix_multiply(4,num_valid_chans,1,(double *)X,(double *)omp,
                    (double *)correction);
    // apply correction to Rx position estimate
    Rx_Pos[0] += correction[0];
    Rx_Pos[1] += correction[1];
    Rx_Pos[2] += correction[2];
    Rx_Clock_Bias = correction[3];
    // if correction is small enough, we're done, exit loop
    tempd = vector_norm(correction);
    if(tempd < 0.001)
    {
      // calculate GDOP
      system_vars.gdop = sqrt(H[0][0] + H[1][1] + H[2][2] + H[3][3]);
      break;
    }
  }  // end of iteration loop

  if (i < max_iterations)
  {
    // solution converged
    // copy final estimate of the receiver position to system variables
    system_vars.recv_pos[0]  = Rx_Pos[0];
    system_vars.recv_pos[1]  = Rx_Pos[1];
    system_vars.recv_pos[2]  = Rx_Pos[2];
    system_vars.recv_pos[3]  = Rx_Clock_Bias;
    // convert to lat/lon/hgt
    wgsxyz2llh(system_vars.recv_pos,system_vars.recv_pos_llh);
    // rotate difference between this pos and reference into NEU
    // use first position as the NEU reference
    if(system_vars.recv_pos_refxyz[0] == 0.0)
    {
      system_vars.recv_pos_refxyz[0] = Rx_Pos[0];
      system_vars.recv_pos_refxyz[1] = Rx_Pos[1];
      system_vars.recv_pos_refxyz[2] = Rx_Pos[2];	
    }
    temp_XYZ[0] = Rx_Pos[0] - system_vars.recv_pos_refxyz[0];
    temp_XYZ[1] = Rx_Pos[1] - system_vars.recv_pos_refxyz[1];
    temp_XYZ[2] = Rx_Pos[2] - system_vars.recv_pos_refxyz[2];
    wgsxyz2neu(temp_XYZ,system_vars.recv_pos_refxyz,system_vars.recv_pos_neu);
    // use this position as the initial guess for the next iteration
    system_vars.initial_pos_guess[0] = Rx_Pos[0];
    system_vars.initial_pos_guess[1] = Rx_Pos[1];
    system_vars.initial_pos_guess[2] = Rx_Pos[2];
  }
  else
  {
    // problem, set return flag
    retval = FAILURE;
  }

  // Set the receiver clock
  if (retval == SUCCESS)
  {
    // steer the receiver clock
    system_vars.recv_time -= Rx_Clock_Bias/NAV_C;
  }

  //####################
  // Velocity Solution #
  //####################

  // Step1, G, Gtrans and H matrix already exists from the 
  // position calculation above.
  if (i < max_iterations)
  {
    // loop through valid measurements
    for (j=0;j<num_valid_chans;j++)
    {
      // calculate Satillite velocity matrix, Step 2
      Vs[j][0] = PVT_Info[j].satvelxyz[0];
      Vs[j][1] = PVT_Info[j].satvelxyz[1];
      Vs[j][2] = PVT_Info[j].satvelxyz[2];
      Vs[j][3] = 0;
      //        pdot[j] = ((NAV_C*dopp_meas[j])/(CARRIER_FREQ)) - PVT_Info[j].satclk_drift*NAV_C;
      pdot[j] = PVT_Info[j].pseudorange_dot;
    }
    matrix_multiply(4, 4, num_valid_chans, (double *)H,
                    (double *)Gtrans, (double *)X);
    for(k=0;k<num_valid_chans;k++)
      tempvX[k] = pdot[k] + G[k][0]*Vs[k][0] + G[k][1]*Vs[k][1] + 
                  G[k][2]*Vs[k][2];
    matrix_multiply(4, num_valid_chans, 1, (double *)X,
                    (double *)tempvX,(double *)velsol);
    system_vars.recv_vel[0] = velsol[0];
    system_vars.recv_vel[1] = velsol[1];
    system_vars.recv_vel[2] = velsol[2];
    system_vars.recv_vel[3] = velsol[3];
  }  // end if max_iterations
  return(retval);
}

/* **************************************************************************/

unsigned int calc_pseudoranges()
{
	unsigned int problem_flags[MAX_CHANNELS], ch_idx;
	unsigned int retval = SUCCESS;
	double temp_time;
	char msg[200];

    // double check we have everything
	for (ch_idx = 0; ch_idx < system_vars.num_channels; ch_idx++)
	{
    // double check that we have a valid subframe and TOW
    problem_flags[ch_idx] = NO;
    if (!c[ch_idx].nav.subframe_start_valid || 
        c[ch_idx].nav.nav_data_state < HAVE_TOW)
      problem_flags[ch_idx] = YES;
	}

  // make channel 0 a reference, for now
  temp_time = c[0].nav.bit_time;

  for (ch_idx = 0; ch_idx < system_vars.num_channels; ch_idx++)
  {
    c[ch_idx].nav.pseudorange_valid = NO;
    if(problem_flags[ch_idx] == NO)
    {
      // after we have set the receiver clock, we can make 
      // pseudorange measurements
      if(system_vars.recv_time_valid)
      {
        // pseudorange taken at some bit boundary
        unsigned bit_codes; // number of complete code transmissions that 
                            // have happened in this nav bit so far
        bit_codes = (c[ch_idx].state_ms - c[ch_idx].track.nav_bit_start) % 
                    MS_PER_NAV_BIT;
        if(bit_codes == 0)
          bit_codes = 20;		
        c[ch_idx].nav.tx_time = c[ch_idx].nav.bit_time + 0.001 * 
                                   (bit_codes + (c[ch_idx].code_prompt - 1) / 
                                   CHIPS_PER_CODE);
        c[ch_idx].nav.pseudorange = (system_vars.recv_time - 
                                     c[ch_idx].nav.tx_time) * NAV_C;
        c[ch_idx].nav.doppler_meas = -1*(c[ch_idx].doppler);
        c[ch_idx].nav.pseudorange_valid = YES;
        c[ch_idx].nav.pseudorange_count++;
        if(ch_idx == 5)
        {
          sprintf(msg, "XXX %.7f %.7f, %d \n", c[ch_idx].nav.bit_time,
                  c[ch_idx].nav.tx_time, bit_codes);
          //XPRINTF(msg);
        }
        c[ch_idx].nav.old_pseudorange = c[ch_idx].nav.pseudorange;
        // solution time is the current bit time
        system_vars.nav_GPS_secs = temp_time;
      }
      else
      {
        c[ch_idx].nav.pseudorange = 0.0;
        c[ch_idx].nav.pseudorange_valid = NO;
      }
//          sprintf(msg,"PRN %d: time=%f, pseudorange=%f\n", c[ch_idx].prn_num, c[ch_idx].nav.pseudorange / NAV_C, c[ch_idx].nav.pseudorange);
//		    XPRINTF(msg);
    }  // if problem_flags
  }  // end channel loop
  system_vars.last_nav_time = c[0].clock;
  return(retval);
}

/* **************************************************************************/

unsigned int WAAS_corrections(unsigned int chan)
{
  unsigned int retval = SUCCESS;
  double range_correction;
  range_correction = CalculateDelayCorrection(chan);
  // apply WAAS correction to pseudoranges
  if(c[chan].nav.pseudorange_valid == YES)
    c[chan].nav.pseudorange += range_correction; 
  return retval;
}

/* ************************************************************** */
//
// Calculates an ionospheric delay correction for a given satellite
// From Section A.4.4.10 of the WAAS MOPS (DO-229D)
// Copywrite 2006, RTCA Inc.
//
/* ************************************************************** */

const double NAV_RADIUS_EARTH = 6378.1363; // km
const double H1               = 350.0;     // km
const double MAX_IGP_DELAY    = 63.876;    // meters

double CalculateDelayCorrection(int chan)
{
  double delay = 0.0;					// delay correction to apply to pseudorange (secs)
  double phi_pp,lambda_pp,phi_pp_deg,lambda_pp_deg;			// lat and lon of IPP
  double psi_pp;						// Earth central angle
  double phi_u,phi_u_deg,lambda_u;	// Rx lat lon
  double SatLLH[3];		// Sat lat lon hgt
  double E,A;				// Sat elevation and azimuth
  double tempd,tempd2,tempd3,tempd4;
  double IGP_points[4][4];  // place for 4 -> [lon,lat,delay,valid_flag]
  int num_IGPs,i;
  double tau_vpp;		// vertical IPP delay
  double x_pp,y_pp;	// 
  double W[4];		// IGP weightings
  double lambda1,lambda2,phi1,phi2;
  double delta_phi_pp,delta_lambda_pp;	
  double F_pp;

  //
  // Determine the Ionospheric Pierce Point IPP
  //

  // convert satellite position to lat/lon
  wgsxyz2llh(c[chan].nav.sat_pos,SatLLH);

  // calculate satellite Az/El with respect to the estimated Rx position, 
  // everything in radians
  phi_u = system_vars.recv_pos_llh[0];
  phi_u_deg = phi_u*R2D;
  lambda_u = system_vars.recv_pos_llh[1];
  CalculateAzEl(system_vars.recv_pos,c[chan].nav.sat_pos,&A,&E);

  // Earth central angle
  tempd = (NAV_RADIUS_EARTH/(NAV_RADIUS_EARTH + H1))*cos(E);
  psi_pp = (M_PI/2.0) - E - asin(tempd);

  // latitude of the IPP
  phi_pp = asin(sin(phi_u)*cos(psi_pp) + cos(phi_u)*sin(psi_pp)*cos(A));

  tempd = tan(psi_pp)*cos(A);
  tempd2 = tan((M_PI/2.0) - phi_u);
  tempd3 = tan(psi_pp)*cos(A + M_PI);
  tempd4 = tan((M_PI/2.0) - phi_u);

  // longitude of the IPP
  if ((phi_u_deg > 70 && tempd > tempd2) || 
      (phi_u_deg < -70 && tempd3 > tempd4))
	  lambda_pp = lambda_u + M_PI - asin((sin(psi_pp)*sin(A))/cos(phi_pp));
  else
    lambda_pp = lambda_u + asin((sin(psi_pp)*sin(A))/cos(phi_pp));

  // Select Ionospheric Grid Points IGP's from File
  phi_pp_deg = phi_pp*R2D;
  lambda_pp_deg = lambda_pp*R2D;
  num_IGPs = ReadWAASFile(phi_pp_deg,lambda_pp_deg,&IGP_points[0][0]);	

  // if we have at least 3 valid points
  if (num_IGPs >= 3 && num_IGPs != FAILURE)
  {
    // sort the IPPs into East-West North-South, a bit ugly
    lambda1 = phi1 = 500.0;
    lambda2 = phi2 = -500.0;
    for (i = 0; i < 4; i++)
    {
      lambda1 = min(lambda1,IGP_points[i][0]);
      lambda2 = max(lambda2,IGP_points[i][0]);
      phi1 = min(phi1,IGP_points[i][1]);
      phi2 = max(phi2,IGP_points[i][1]);
    }

    delta_phi_pp = phi_pp_deg - phi1;
    delta_lambda_pp = lambda_pp_deg - lambda1;
    x_pp = delta_lambda_pp/(lambda2 - lambda1); 
    y_pp = delta_phi_pp/(phi2 - phi1);

    // Interpolate the grid points to the IPP
    if (num_IGPs == 4)
    {
      // IPP in the middle of 4 valid grid points
      // calculate weights, note: only works if IPP is withing 85deg lat
      W[0] = x_pp*y_pp;
      W[1] = (1 - x_pp)*y_pp;
      W[2] = (1 - x_pp)*(1 - y_pp);
      W[3] = x_pp*(1 - y_pp);
      // calculate vertical IPP delay
      tau_vpp = 0.0;
      for(i = 0; i < 4; i++)
        tau_vpp +=  W[i]*IGP_points[i][2];  // weight * delay
    }
    else if(num_IGPs == 3)
    {
      // IPP in the middle of 3 valid grid points
      W[0] = y_pp;
      W[1] = 1 - x_pp - y_pp;
      W[2] = x_pp;
      // calculate vertical IPP delay
      tau_vpp = 0.0;
      for(i = 0; i < 3; i++)
        tau_vpp +=  W[i]*IGP_points[i][2];  // weight * delay
    }
    else
      return delay;  // 0.0

    // Compute Slant Delay
    // compute obliquity factor
    tempd = (NAV_RADIUS_EARTH*cos(E))/(NAV_RADIUS_EARTH + H1);
    tempd2 = 1 - tempd*tempd;
    F_pp = 1.0/sqrt(tempd2); 
    // compute pseudorange connection
    delay = -1*F_pp*tau_vpp;
    // whew ...
  }
  return delay;
}

/* **************************************************************
 * Reads upto 4 IGP points from the specified external file
 **************************************************************** */

int ReadWAASFile(double phi_pp_deg,double lambda_pp_deg,double *IGP_points)
{
  double temp_lon,temp_lat,temp_delay;
  double tempd;
  int point_count = 0;
  int corners_found = 0;
//
// Currently, this routine only looks for the points within a 5 degree grid
// around the IPP. so... not all of the possibities from Appendix A are covered.
// IGPs at very high latitudes may be problems
// on va voir
//
  if(system_vars.waas_igp_file == NULL)
    system_vars.waas_igp_file = fopen(system_vars.WAASfilename,"r");
  if(!system_vars.waas_igp_file) 
    return FAILURE;
  else
  {
	  rewind(system_vars.waas_igp_file);
    while(!feof(system_vars.waas_igp_file))
    {
      // read in longitude, deg
      VERIFY_IO(fscanf(system_vars.waas_igp_file," %lf",&temp_lon), 1);
      // read in latitude, deg
      VERIFY_IO(fscanf(system_vars.waas_igp_file," %lf",&temp_lat), 1);
      // read in delay, meters
      VERIFY_IO(fscanf(system_vars.waas_igp_file," %lf",&temp_delay), 1);
      // is this point within 5 deg of IPP longitude (lambda)
      tempd = fabs(lambda_pp_deg - temp_lon); 
      if(tempd < 5.0)
      {
        // is this point within 5 deg of IPP latitude (phi)
        tempd = fabs(phi_pp_deg - temp_lat); 
        if(tempd < 5.0)
        {
          // this is a corner
          corners_found++;
          // is this a valid delay
          if(temp_delay <= MAX_IGP_DELAY)
          {
            // save the point
            IGP_points[point_count*4 + 0] = temp_lon;		// lon
            IGP_points[point_count*4 + 1] = temp_lat;		// lat
            IGP_points[point_count*4 + 2] = temp_delay;	// delay
            IGP_points[point_count*4 +3] = 1;		// valid
            // increment to next point
            point_count++;
          } 
        }  // end lat check
      }  // end lon check
      // when we have found all 4 grid points surounding the IPP
      if(corners_found >= 4)
        return point_count;
    }  // end while
  } // end else
  return point_count;
}

/***************************************************************************/

void CalculateAzEl(double RxXYZ[3], double SatXYZ[3], double *A,double *E)
{
  double tempd,tempd2;
  double temp_XYZ[3];
  double NEU[3];
  // local vector to rotate, from Rx to Sat
  temp_XYZ[0] = SatXYZ[0] - RxXYZ[0];
  temp_XYZ[1] = SatXYZ[1] - RxXYZ[1];
  temp_XYZ[2] = SatXYZ[2] - RxXYZ[2];
  // rotate this vector into NEU
  wgsxyz2neu(temp_XYZ,RxXYZ,NEU);
  // compute elevation in local NEU frame
  tempd = sqrt(NEU[0]*NEU[0] + NEU[1]*NEU[1]);
  tempd2 = NEU[2]/tempd;
  *E = atan(tempd2);
  // compute azimuth in local NEU frame
  *A = atan2(NEU[1],NEU[0]);
}

/***************************************************************************/

void nav_init()
{
  unsigned int ch_idx;
  for (ch_idx = 0; ch_idx < system_vars.num_channels; ch_idx++)
  {
    int sf_idx;
    struct channel *ch = &c[ch_idx];
    ch->nav.subframe_write_pos = 0;
    for (sf_idx = 0; sf_idx < SUBFRAME_LENGTH; sf_idx++)
      ch->nav.subframe[sf_idx] = 0;
    ch->nav.subframe_state = SF_SEARCHING;
    ch->nav.num_preamb_cand = 0;
    ch->nav.nav_data_state = 0;
    for (sf_idx = 0; sf_idx < 6; sf_idx++)
      ch->nav.subframe_valid[sf_idx] = 0;
    ch->nav.subframe_start_valid = 0;
    ch->nav.subframe_start_clock = 0;
    ch->nav.recv_pos[0] = ch->nav.recv_pos[1] = 0;
    ch->nav.recv_pos[2] = ch->nav.recv_pos[3] = 0;
    ch->nav.first_subframe_flag = 0;
  }
  system_vars.recv_pos_refxyz[0] = 0.0;
  system_vars.recv_pos_refxyz[1] = 0.0;
  system_vars.recv_pos_refxyz[2] = 0.0;
}

/****************************************************************************
* INT16U GetSVInfo(s_SV_Info *sv)
*
* Interface to intrpsp3.cpp functions to determine satellite pos/vel and 
* clock bias from sp3 files.
*
* Dec 2007, STG
*
* Inputs:       sv->prn, sv->week, sv->TOW
* Outputs:      sv->posxyz, sv->velxyz, sv->clk_bias
*
****************************************************************************/
int GetSVInfo(struct s_SV_Info *sv, char *infileXname)
{
  int retval = SUCCESS;
  //int prnNum;
  char prnNum[10];
  string prnNum_str;
  DateTime currEpoch;
  double PosVel[10];
  string orbfile;
  GPSTime x;
  SP3cFile mysp3;

  // this hack prevents the occasional crash on array copying (on some platforms)
  sprintf(prnNum,"G%2d",sv->prn);
  if(sv->prn < 10)
	prnNum[1] = '0';  // I hate strings
  for(int i = 0; i < 3; i++)
    prnNum_str += prnNum[i];

  x.GPSWeek = sv->week;
  x.secsOfWeek = sv->TOW;
  currEpoch = DateTime(x);
  mysp3.setPathFilenameMode(infileXname);
  mysp3.readHeader();
  retval = (int) mysp3.getSVPosVel( currEpoch, prnNum_str, PosVel);
  if(retval == 0)
  {
    sv->posxyz[0] = PosVel[0]*1000;  // m
    sv->posxyz[1] = PosVel[1]*1000;
    sv->posxyz[2] = PosVel[2]*1000;
    sv->clk_bias  = PosVel[3]/1e6;   // sec
    sv->velxyz[0] = PosVel[4]*1000;  // m/s
    sv->velxyz[1] = PosVel[5]*1000;
    sv->velxyz[2] = PosVel[6]*1000;
    sv->clk_drift = PosVel[7]/1e6;   // sec/sec
  }
  else
  {
    sv->posxyz[0] = 0;
    sv->posxyz[1] = 0;
    sv->posxyz[2] = 0;
    sv->clk_bias  = 0;
    sv->velxyz[0] = 0;
    sv->velxyz[1] = 0;
    sv->velxyz[2] = 0;
    sv->clk_drift = 0;
  }
  return retval;
}

/****************************************************************************
* Calculate satellite position, velocity and clock offset from ephemeris
****************************************************************************/
unsigned int ProcessEphemeris(unsigned int week, double secs, 
                              unsigned int sv, nav_info_t *Eph)
{
  unsigned int retval = SUCCESS;
  double tempd1 = 0, tempd2, tempd3;
  double tdiff;
  double a;               // semi major axis
  double ma,ma_dot;       // mean anomoly and first derivative
  double ea,ea_dot,ea_old;   // eccentric anomoly, first deriv, iteration var
  double einstein;        // relativistic correction
  double al,al_dot;       // argument of lattitude and first derivative
  double cal,cal_dot;     // corrected argument of lattitude and first deriv
  double r,r_dot;         // radius and first derivative
  double inc,inc_dot;     // inclination and first derivative
  double x,x_dot,y,y_dot; // position in orbital plan and first derivatives
  double om,om_dot;     // omega and first derivatives

  // seconds from the toe
  tdiff = secs - (double) Eph->toe;

  if (tdiff > 302400.0)
    tdiff -= 604800.0;
  else if(tdiff < -302400.0)
    tdiff += 604800.0;

  a = Eph->sqrta*Eph->sqrta;
  ma_dot = sqrt(NAV_GM/(a*a*a)) + Eph->dn;
  ma = Eph->m0 + ma_dot*tdiff;

  ea = ma;
  ea_old = ea + 1;

  while(fabs(ea-ea_old) > 1.0E-14)
  {
    ea_old = ea;
    tempd1 = 1.0 - Eph->ecc*cos(ea_old);
    ea = ea + (ma-ea_old+Eph->ecc*sin(ea_old))/tempd1;
  }
  ea_dot = ma_dot/tempd1;

  einstein = -4.442807633E-10*Eph->ecc*Eph->sqrta*sin(ea);
  tempd2 = sqrt(1.0 - Eph->ecc*Eph->ecc);

  al = atan2(tempd2*sin(ea),cos(ea)-Eph->ecc) + Eph->w;
  al_dot = tempd2*ea_dot/tempd1;

  cal = al + Eph->cus*sin(2.0*al) + Eph->cuc*cos(2.0*al);
  cal_dot = al_dot*(1.0 + 2.0*(Eph->cus*cos(2.0*al) - Eph->cuc*sin(2.0*al)));

  r = a*tempd1 + Eph->crc*cos(2.0*al) + Eph->crs*sin(2.0*al);
  r_dot = a*Eph->ecc*sin(ea)*ea_dot + 2.0*al_dot*(Eph->crs*cos(2.0*al) - 
                                                  Eph->crc*sin(2.0*al));

  inc = Eph->inc + Eph->inc_dot*tdiff + Eph->cic*cos(2.0*al) + 
        Eph->cis*sin(2.0*al);
  inc_dot = Eph->inc_dot + 2.0*al_dot*(Eph->cis*cos(2.0*al) - 
                                       Eph->cic*sin(2.0*al));

  x = r*cos(cal);
  y = r*sin(cal);
  x_dot = r_dot*cos(cal) - y*cal_dot;
  y_dot = r_dot*sin(cal) + x*cal_dot;

  om_dot = Eph->omegadot - NAV_OMEGAE_DOT;
  om = Eph->omega0 + tdiff*om_dot - NAV_OMEGAE_DOT*Eph->toe;

  /* Compute the satellite's position. */
  Eph->sat_pos[0] = x*cos(om) - y*cos(inc)*sin(om);
  Eph->sat_pos[1] = x*sin(om) + y*cos(inc)*cos(om);
  Eph->sat_pos[2] = y*sin(inc);

  tempd3 = y_dot*cos(inc) - y*sin(inc)*inc_dot;

  /* Compute the satellite's velocity. */
  Eph->sat_vel[0] = -om_dot*Eph->sat_pos[1] + x_dot*cos(om) - tempd3*sin(om);
  Eph->sat_vel[1] = om_dot*Eph->sat_pos[0] + x_dot*sin(om) + tempd3*cos(om);
  Eph->sat_vel[2] = y*cos(inc)*inc_dot + y_dot*sin(inc);

  // seconds from the toc
  tdiff = secs - Eph->toc;

  if(tdiff > 302400.0)
    tdiff -= 604800.0;
  else if(tdiff < -302400.0)
    tdiff += 604800.0;

  // satellite clock terms
  Eph->clock_corr = Eph->af0 + tdiff*(Eph->af1 + tdiff*Eph->af2) + 
                    einstein - Eph->tgd;
  Eph->clock_drift = Eph->af1 + 2.0*tdiff*Eph->af2;

  return retval;
}

