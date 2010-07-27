/*
* Copyright (c) 2010, Scott Gleason and Kim Chang
* All rights reserved.
*
* These routines will estimate the receiver position using 1ms code phase measurements
* It requires a decent estimate of the rx postion (<150km) and time (<60s).

* Details can be found in the following publications,
*
* Peterson, B., R. Hartnett, and G. Ottman, GPS Receiver Structures for the Urban Canyon, Proc. of the ION-GPS, 1995.
* Lannelongue, S.; Pablos, P. Fast Acquisition Techniques For GPS Receivers, Proc. of the ION-GPS, 1998.
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

extern void invert(unsigned actualsize,double *data);
extern void invert4x4(double A[4][4], double Ainv[4][4]);
extern double vector_norm(double v[3]);
extern void vector_subtract(double a[3], double b[3],double c[3]);
extern double  vector_dot_product (double *a, double *b);
extern void matrix_multiply(unsigned rowsA,unsigned columnsA,unsigned columnsB,double *A,double *B,double *C);
extern void matrix_transpose(unsigned rows,unsigned columns,double *A,double *B);
extern void wgsllh2xyz(double *llh,double *xyz);
extern void wgsxyz2llh(double *xyz,double *llh);
extern double predicted_doppler(int gpsweek, double gpssecs,int prn,double *rxpos,double *rxvel);


FILE* debug_file;
FILE* nav_file2;
s_misc_info misc_info;

//* *******************************************************************************
//  SingleEpochPosition
//* *******************************************************************************

int SingleEpochPosition(void){
int retval = 0;
unsigned ch_idx;
s_nav_packet nav_pkt;

          memset(&nav_pkt,0,sizeof(struct s_nav_packet));
          nav_pkt.gps_week = system_vars.estimate_gpsweek;
          nav_pkt.gps_seconds = system_vars.estimate_gpssecs;

          // don't need this but ... if its in the file
          nav_pkt.pos_truth[0] = system_vars.estimate_wgs84_pos[0];
          nav_pkt.pos_truth[1] = system_vars.estimate_wgs84_pos[1];
          nav_pkt.pos_truth[2] = system_vars.estimate_wgs84_pos[2];

          nav_pkt.time_init = system_vars.estimate_gpssecs;
          if(system_vars.position_status == HAVE_WGS84_FILE_ESTIMATE){
            nav_pkt.pos_init[0] = system_vars.estimate_wgs84_pos[0];
            nav_pkt.pos_init[1] = system_vars.estimate_wgs84_pos[1];
            nav_pkt.pos_init[2] = system_vars.estimate_wgs84_pos[2];
          }else if(system_vars.position_status == HAVE_DOPPLER_POS_ESTIMATE){
            nav_pkt.pos_init[0] = system_vars.recv_dopp_pos[0];
            nav_pkt.pos_init[1] = system_vars.recv_dopp_pos[1];
            nav_pkt.pos_init[2] = system_vars.recv_dopp_pos[2];
          }else{
            // HAVE_SNAPSHOT_POS_ESTIMATE
            nav_pkt.pos_init[0] = system_vars.recv_snapshot_pos[0];
            nav_pkt.pos_init[1] = system_vars.recv_snapshot_pos[1];
            nav_pkt.pos_init[2] = system_vars.recv_snapshot_pos[2];
            nav_pkt.time_init = system_vars.recv_snapshot_time + system_vars.PVT_INTERVAL;
          }

          nav_pkt.num_meas = system_vars.num_channels;

          for (ch_idx = 0; ch_idx < nav_pkt.num_meas; ch_idx++)
          {
              nav_pkt.Meas[ch_idx].prn = c[ch_idx].prn_num;
              nav_pkt.Meas[ch_idx].ms_range = (c[ch_idx].code_prompt-1)/1023.0;
              nav_pkt.Meas[ch_idx].doppler = -1*(c[ch_idx].doppler);  // -1*
          }

          retval = CalcSingleEpochPosition(&nav_pkt,1);

return(retval);
}

//* *************************************************
//  Calculate SingleEpochPosition
// *************************************************

int CalcSingleEpochPosition(s_nav_packet *pkt,int signal){
int retval = 0;
s_PVT_Info2  PVT_Info[MAX_CHANNELS];
double Rx_Pos[3];
double Rx_Vel[3];
double Rx_Clock_Bias;
unsigned int max_iterations,i,j,k;
double omp[MAX_CHANNELS];
double G[MAX_CHANNELS][5];
double Gtrans[5][MAX_CHANNELS];
double GtG[5][5];
double v0vk_mag,tau,wEtau;
double rotm[3][3];
double vk_new[3];
double tempv[3];
double los[3];
double correction[5];
double H[5][5];
double X[5][MAX_CHANNELS];
unsigned num_valid_chans = pkt->num_meas;
unsigned num_valid_meas = 0;
double Rx_Pos_Truth[3];
double tg_Truth;
unsigned int week;
double tg;  // coarse GPS time, unknown
double ri[MAX_CHANNELS],si[MAX_CHANNELS],pi[MAX_CHANNELS],p_dot[MAX_CHANNELS];
double sat_clk_m[MAX_CHANNELS],sat_clk_ms[MAX_CHANNELS];
double ms_range[MAX_CHANNELS];
double tx_time[MAX_CHANNELS],pi2[MAX_CHANNELS];
double tempd,tempd2;
double temptime;
double tempbits[MAX_CHANNELS];
double sat_clks_ms[MAX_CHANNELS];
s_SV_Info sv_info;
double temp_XYZ[3];

Rx_Pos_Truth[0] = pkt->pos_truth[0];
Rx_Pos_Truth[1] = pkt->pos_truth[1];
Rx_Pos_Truth[2] = pkt->pos_truth[2];
tg_Truth = pkt->gps_seconds;
week = pkt->gps_week;

Rx_Pos[0] = pkt->pos_init[0];   // WGS84 X
Rx_Pos[1] = pkt->pos_init[1];   // WGS84 Y
Rx_Pos[2] = pkt->pos_init[2];   // WGS84 Z
tg = pkt->time_init;

Rx_Vel[0] = 0.0;
Rx_Vel[1] = 0.0;
Rx_Vel[2] = 0.0;

vector_subtract(Rx_Pos,Rx_Pos_Truth,tempv);
misc_info.pos_diff_mag = vector_norm(tempv);
misc_info.time_diff_mag = tg - tg_Truth;
misc_info.final_pos_diff_mag = -99;
misc_info.final_time_diff_mag = -88;
misc_info.final_iterations = 101;

// extract valid measurements at this time epoch
for (i=0;i<num_valid_chans;i++){


		PVT_Info[num_valid_meas].prn = pkt->Meas[i].prn;
		PVT_Info[num_valid_meas].TOA_m = pkt->Meas[i].ms_range; // modulo 1ms
        PVT_Info[num_valid_meas].doppler = pkt->Meas[i].doppler;

        // calculate satellite positions at this time
        sv_info.prn = PVT_Info[num_valid_meas].prn;
        sv_info.week = week;
        sv_info.TOW = tg;
        retval = GetSVInfo(&sv_info,system_vars.IGSfilename);

        if(retval == 0){
                num_valid_meas++;
        }
}

if(num_valid_meas < 5)
    return(retval);

//num_valid_meas = 6;

//####################
// Position Solution #
//####################

// **************************************************************
// start of iteration loop
tempd = 0;
max_iterations = 10;
//max_iterations = 1;
misc_info.tg_correction = 1e10;
for (i=0;i<max_iterations;i++){

	// calculate predicted ms ranges for all measurements
    for (j=0;j<num_valid_meas;j++){

        // calculate satellite positions at tg
        sv_info.prn = PVT_Info[j].prn;
        sv_info.week = week;
        sv_info.TOW = tg;
        retval = GetSVInfo(&sv_info,system_vars.IGSfilename);

        // calculate estimated travel time using tg
        vector_subtract(Rx_Pos,sv_info.posxyz,tempv);
        tempd = vector_norm(tempv);

        // recalculate satellite postions at tg - r/c
        tempd2 = tg - (tempd/SPEED_OF_LIGHT);
        sv_info.TOW = tempd2;
        retval = GetSVInfo(&sv_info,system_vars.IGSfilename);

        if(retval == 0){
                PVT_Info[j].satposxyz[0] = sv_info.posxyz[0];
                PVT_Info[j].satposxyz[1] = sv_info.posxyz[1];
                PVT_Info[j].satposxyz[2] = sv_info.posxyz[2];
                PVT_Info[j].satclk_bias  = sv_info.clk_bias;
                PVT_Info[j].satvelxyz[0] = sv_info.velxyz[0];
                PVT_Info[j].satvelxyz[1] = sv_info.velxyz[0];
                PVT_Info[j].satvelxyz[2] = sv_info.velxyz[0];
                PVT_Info[j].satclk_drift = sv_info.clk_drift;
        }else{
            return(retval);
        }

        // The satellite positions need to be corrected for Earth Rotation during the transmission time.
        vector_subtract(Rx_Pos,PVT_Info[j].satposxyz,tempv);
        v0vk_mag = vector_norm(tempv);
        tau = v0vk_mag/SPEED_OF_LIGHT;
		wEtau = OMEGAE_DOT*tau;

		rotm[0][0] = cos(wEtau);
		rotm[0][1] = sin(wEtau);
		rotm[0][2] = 0.0;
		rotm[1][0] = -sin(wEtau);
		rotm[1][1] = cos(wEtau);
		rotm[1][2] = 0.0;
		rotm[2][0] = 0.0;
		rotm[2][1] = 0.0;
		rotm[2][2] = 1.0;

        matrix_multiply(3,3,1,(double *)rotm,(double *)PVT_Info[j].satposxyz,(double *)vk_new);  // result in vk_new

        // calculate predicted pseudorange for this satellite
        vector_subtract(PVT_Info[j].satposxyz,Rx_Pos,tempv);
        vector_subtract(Rx_Pos,vk_new,tempv);
        ri[j] = vector_norm(tempv);

        // predicted pseudorange in milliseconds
        ms_range[j] = (ri[j]/(0.001*SPEED_OF_LIGHT));

        // code phase based, millisecond remainder
        si[j] = PVT_Info[j].TOA_m;

        // pseudorange rate
        p_dot[j] = PVT_Info[j].doppler*L1_WAVELENGTH;

		// line of sight unit vector
		los[0] = tempv[0]/ri[j];
		los[1] = tempv[1]/ri[j];
		los[2] = tempv[2]/ri[j];

		// build geometry matrix
		G[j][0] = 1*los[0];
		G[j][1] = 1*los[1];
		G[j][2] = 1*los[2];
        G[j][3] = SPEED_OF_LIGHT;
		G[j][4] = 1*p_dot[j];

    }


    // loop through channels, calculate pseudoranges
    for (j=0;j<num_valid_meas;j++){

            // L1 1ms ranges
            pi[j] = (floor(ms_range[j])*1.0e-3 + si[j]*1.0e-3)*SPEED_OF_LIGHT;
            sat_clk_ms[j] = PVT_Info[j].satclk_bias*1000;
            sat_clk_m[j] = PVT_Info[j].satclk_bias*SPEED_OF_LIGHT;
            pi[j] += PVT_Info[j].satclk_bias*SPEED_OF_LIGHT;

            sat_clks_ms[j] = PVT_Info[j].satclk_bias*1000;

            tempbits[j] = round(100.0 - ms_range[j] + (si[0] - si[j]) + sat_clks_ms[j]);

            temptime = tg - 5.0;

            tx_time[j] = (temptime + tempbits[j]*1.0e-3 + si[j]*1.0e-3);
            pi2[j] = (tg - tx_time[j])*SPEED_OF_LIGHT;
            pi2[j] += PVT_Info[j].satclk_bias*SPEED_OF_LIGHT;

            // observed minus predicted PR
            omp[j] = (pi2[j] - ri[j]);

            if(omp[j] > 149000){
                omp[j] -= 1.0e-3*SPEED_OF_LIGHT;
            }

            if(omp[j] < -149000){
                omp[j] += 1.0e-3*SPEED_OF_LIGHT;
            }

    }  // end j

	// solve for corrections
    matrix_transpose(num_valid_meas,5,(double *)G,(double *)Gtrans);
    matrix_multiply(5,num_valid_meas,5,(double *)Gtrans,(double *)G,(double *)GtG);
    invert(5,(double *) GtG);  // result in GtG
    memcpy(H,GtG,sizeof(GtG));  // copy to H
    matrix_multiply(5,5,num_valid_meas,(double *)H,(double *)Gtrans,(double *)X);
    matrix_multiply(5,num_valid_meas,1,(double *)X,(double *)omp,(double *)correction);

	// apply correction to Rx position estimate
    Rx_Pos[0] += correction[0];
    Rx_Pos[1] += correction[1];
    Rx_Pos[2] += correction[2];
    Rx_Clock_Bias = correction[3];

    tg += correction[4];

    misc_info.tg_correction = correction[4];

    misc_info.pos_correction_mag = vector_norm(correction);

	// if correction is small enough, we're done, exit loop
    if((misc_info.pos_correction_mag < POS_TOL) && (misc_info.tg_correction < TG_TOL)){

        // calculate PDOP
        misc_info.pdop = sqrt(H[0][0] + H[1][1] + H[2][2]);

        break;

    }

}

// end of iteration loop
// **************************************************************

for (k=0;k<num_valid_meas;k++){misc_info.debug[k] = omp[k];}
update_debug_log(num_valid_meas);

if (i < max_iterations){

    // success, save result
    retval = SUCCESS;

    misc_info.recv_pos[0]  = Rx_Pos[0];
    misc_info.recv_pos[1]  = Rx_Pos[1];
    misc_info.recv_pos[2]  = Rx_Pos[2];
    misc_info.clock_bias  = Rx_Clock_Bias;
    misc_info.tg  = tg;

    system_vars.recv_snapshot_pos[0] = Rx_Pos[0];
    system_vars.recv_snapshot_pos[1] = Rx_Pos[1];
    system_vars.recv_snapshot_pos[2] = Rx_Pos[2];
    system_vars.recv_snapshot_pos[3] = Rx_Clock_Bias;
    system_vars.recv_snapshot_time = tg;

    system_vars.recv_pos[0]  = Rx_Pos[0];
    system_vars.recv_pos[1]  = Rx_Pos[1];
    system_vars.recv_pos[2]  = Rx_Pos[2];
    system_vars.recv_pos[3]  = Rx_Clock_Bias;
    system_vars.recv_time = tg;

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

    system_vars.position_status = HAVE_SNAPSHOT_POS_ESTIMATE;

}else{

    misc_info.recv_pos[0]  = 99.0;
    misc_info.recv_pos[1]  = 99.0;
    misc_info.recv_pos[2]  = 99.0;
    misc_info.clock_bias  = 0.0;
    misc_info.tg  = 0.0;

    // problem, set return flag
    retval = PROBLEM;

}

// see how close to truth we got
vector_subtract(Rx_Pos,Rx_Pos_Truth,tempv);
misc_info.final_pos_diff_mag = vector_norm(tempv);
misc_info.final_time_diff_mag = tg - tg_Truth;
misc_info.final_iterations = i;

//update_nav_log2();

return(retval);
}

//****************************************************************************
// INT16U GetSVInfo(s_SV_Info *sv)
//
// Interface to intrpsp3.cpp functions to determine satellite pos/vel and
// clock bias from sp3 files.
//
// Dec 2007, STG
//
// Inputs:       prn, gpsweek, gpssecs, filename
// Outputs:      pos X,Y,Z, clkbias, Vel X,Y,Z clk drift
//
//****************************************************************************
int GetSVInfo2(double *PosVelClk, int prn, int gpsweek, double gpssecs, char *infileXname)
{
  int retval = 0;
  //int prnNum;
  char prnNum[10];
  string prnNum_str;
  DateTime currEpoch;
  double PosVel[10];
  string orbfile;
  GPSTime x;
  SP3cFile mysp3;

  // this hack prevents the occasional crash on array copying (on some platforms)
  sprintf(prnNum,"G%2d",prn);
  if(prn < 10)
	prnNum[1] = '0';  // I hate strings
  for(int i = 0; i < 3; i++)
    prnNum_str += prnNum[i];

  x.GPSWeek = gpsweek;
  x.secsOfWeek = gpssecs;
  currEpoch = DateTime(x);
  mysp3.setPathFilenameMode(infileXname);
  mysp3.readHeader();
  retval = (int) mysp3.getSVPosVel( currEpoch, prnNum_str, PosVel);
  if(retval == 0)
  {
    PosVelClk[0] = PosVel[0]*1000;  // m
    PosVelClk[1] = PosVel[1]*1000;
    PosVelClk[2] = PosVel[2]*1000;
    PosVelClk[3] = PosVel[3]/1e6;   // sec
    PosVelClk[4] = PosVel[4]*1000;  // m/s
    PosVelClk[5] = PosVel[5]*1000;
    PosVelClk[6] = PosVel[6]*1000;
    PosVelClk[7] = PosVel[7]/1e6;   // sec/sec
  }
  else
  {
    PosVelClk[0] = 0;
    PosVelClk[1] = 0;
    PosVelClk[2] = 0;
    PosVelClk[3] = 0;
    PosVelClk[4] = 0;
    PosVelClk[5] = 0;
    PosVelClk[6] = 0;
    PosVelClk[7] = 0;
  }
  return retval;
}

// **************************************************************
//  update_debug_log
// **************************************************************
void update_debug_log(int num_points)
{
int i;

    if(debug_file == NULL){
        debug_file = fopen("debug_log2.dat","w");
    }

    if(debug_file != NULL){

      if(num_points < 60){

        for(i=0;i<num_points;i++){
        	fprintf(debug_file, "%15.9f ",misc_info.debug[i]);
        }
        fprintf(debug_file, "\n");

      }else{
        fprintf(debug_file, "XXX\n");
      }

    }
}


// **************************************************************
//  update_nav_log2
// **************************************************************
void update_nav_log2(void)
{

    if(system_vars.nav_log2 == NULL){
        system_vars.nav_log2 = fopen("nav_log2.dat","w");
    }

    if(system_vars.nav_log2 != NULL){

        fprintf(system_vars.nav_log2, "%15.9f ",misc_info.tg);
        fprintf(system_vars.nav_log2, "%15.9f %15.9f %15.9f ",misc_info.recv_pos[0],misc_info.recv_pos[1],misc_info.recv_pos[2]);
        fprintf(system_vars.nav_log2, "%.5f ",misc_info.clock_bias);
        fprintf(system_vars.nav_log2, "%.5f ",misc_info.final_pos_diff_mag);
        fprintf(system_vars.nav_log2, "%.5f ",misc_info.final_time_diff_mag);
        fprintf(system_vars.nav_log2, "%d ",misc_info.final_iterations);
        fprintf(system_vars.nav_log2, "\n");

    }
}


/****************************************************************************
* Function: double predicted_doppler(double time,double *rxpos)
****************************************************************************/
double predicted_doppler(int gpsweek, double gpssecs,int prn,double *rxpos,double *rxvel)
{
double retval;
double sv_info_IGS[10];
double a,b;
double tempd1;
double v[3];
double unit1[3];
double satvel[3];
double dopp;

    // calculate satellite positions at this time
    retval = GetSVInfo2(sv_info_IGS,prn,gpsweek,gpssecs,(char*)"data/igu14794_06.sp3");

    /* calculate Doppler offset measurement */

    // determine unit vectors from Rx to SV
    vector_subtract(sv_info_IGS,rxpos,v);

    tempd1 = vector_norm(v);
    unit1[0] = v[0]/tempd1;
    unit1[1] = v[1]/tempd1;
    unit1[2] = v[2]/tempd1;

    // project velocities onto this unit vector, wrt Rx
    satvel[0] = -1*sv_info_IGS[4];
    satvel[1] = -1*sv_info_IGS[5];
    satvel[2] = -1*sv_info_IGS[6];
    a = vector_dot_product(rxvel,unit1);
    b = vector_dot_product(satvel,unit1);

    // calculate Doppler measurement, we are ignoring doppler due to sat clock ...
	dopp = -1*((a + b)*CARRIER_FREQ)/SPEED_OF_LIGHT;

return(dopp);
}


