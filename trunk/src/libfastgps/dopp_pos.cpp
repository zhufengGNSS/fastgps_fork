/*
* Copyright (c) 2010, Scott Gleason and Alex Briere
* All rights reserved.
*
* These routines will estimate the receiver position using Doppler frequency measurements
* It requires a decent estimate of the time (<60s). The position estimated is accurate enough
* to initialise the snapshot position routines, but probably not much else.

* Details can be found in the following publication,
*
* Jonathan Hill, The principle of snapthot navigation solution based on Doppler shift, ION GPS 2001 Sep.2001.
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

extern int GetSVInfo2(double *PosVelClk, int prn, int gpsweek, double gpssecs, char *infileXname);

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

FILE* debug_file3;
s_misc_info_dopp misc_info_dopp;

//* *******************************************************************************
//  Doppler Position
//* *******************************************************************************

int DopplerPosition(void){
int retval = 0;
s_nav_packet_dopp nav_pkt_dopp;
unsigned ch_idx;

    memset(&nav_pkt_dopp,0,sizeof(struct s_nav_packet_dopp));

    nav_pkt_dopp.gps_week = system_vars.estimate_gpsweek;
    nav_pkt_dopp.gps_seconds = system_vars.estimate_gpssecs;
    nav_pkt_dopp.pos_truth[0] = system_vars.estimate_wgs84_pos[0];
    nav_pkt_dopp.pos_truth[1] = system_vars.estimate_wgs84_pos[1];
    nav_pkt_dopp.pos_truth[2] = system_vars.estimate_wgs84_pos[2];
    nav_pkt_dopp.num_meas = system_vars.num_channels; // ???

    for (ch_idx = 0; ch_idx < nav_pkt_dopp.num_meas; ch_idx++)
    {
        nav_pkt_dopp.Meas[ch_idx].prn = c[ch_idx].prn_num;
        nav_pkt_dopp.Meas[ch_idx].doppler = (c[ch_idx].doppler);  // no -1?
    }

    retval = CalcDopplerPosition(&nav_pkt_dopp);

return(retval);
}

//* *******************************************************************************
//  Doppler Position Solution
//  Ref: "The Principal of a Snapshot Navigation Solution Based on Doppler Shift"
//	 J. Hill, Proceedings of the ION GPS 2001
//* *******************************************************************************

int CalcDopplerPosition(s_nav_packet_dopp *pkt){
int retval = 0;
double Xc[4],Xc_last[4];
unsigned num_valid_chans = pkt->num_meas;
unsigned num_valid_meas = 0;
s_PVT_Info_dopp  PVT_Info[MAX_CHANNELS];
s_SV_Info sv_info;
unsigned int max_iterations,i,j;
double gpstime_sec;
unsigned int gpstime_week;
double correction,Ym,tempd,err_last;
double err = 0;
double A[MAX_CHANNELS][4];
double tempv[3];
double Rx_Pos_Truth[3];
double B[MAX_CHANNELS];
double AtB[MAX_CHANNELS];
double W[MAX_CHANNELS];
double AtA[4][4];
double invAtA[4][4];
double Atrans[4][MAX_CHANNELS];

Rx_Pos_Truth[0] = pkt->pos_truth[0];
Rx_Pos_Truth[1] = pkt->pos_truth[1];
Rx_Pos_Truth[2] = pkt->pos_truth[2];

gpstime_sec = pkt->gps_seconds;
gpstime_week = pkt->gps_week;

// Rx position and clock bias initial guess
Xc[0] = 0;Xc[1] = 0;Xc[2] = 0;Xc[3] = 0;
Xc_last[0] = 0;Xc_last[1] = 0;Xc_last[2] = 0;Xc_last[3] = 0;

// extract valid measurements at this time epoch
for (i=0;i<num_valid_chans;i++){


	PVT_Info[num_valid_meas].prn = pkt->Meas[i].prn;
    PVT_Info[num_valid_meas].doppler = pkt->Meas[i].doppler;

    W[num_valid_meas] = (PVT_Info[num_valid_meas].doppler*SPEED_OF_LIGHT)/CARRIER_FREQ;  // no -1, speed from sat to Rx

    // calculate satellite positions at this time
    sv_info.prn = PVT_Info[num_valid_meas].prn;
    sv_info.week = gpstime_week;
    sv_info.TOW = gpstime_sec;
    retval = GetSVInfo(&sv_info,system_vars.IGSfilename);

    if(retval == 0){
        num_valid_meas++;
    }
}

if(num_valid_meas < 5)
    return(retval);


//####################
// Position Solution #
//####################

// calculate A matrix at gpstime_sec
for (j=0;j<num_valid_meas;j++){

    // calculate satellite velocities at best gpstime_sec
    sv_info.prn = PVT_Info[j].prn;
    sv_info.week = gpstime_week;
    sv_info.TOW = gpstime_sec;
    retval = GetSVInfo(&sv_info,system_vars.IGSfilename);

    if(retval == 0){
        PVT_Info[j].satposxyz[0] = sv_info.posxyz[0];
        PVT_Info[j].satposxyz[1] = sv_info.posxyz[1];
        PVT_Info[j].satposxyz[2] = sv_info.posxyz[2];
        PVT_Info[j].satclk_bias  = sv_info.clk_bias;
        PVT_Info[j].satvelxyz[0] = sv_info.velxyz[0];
        PVT_Info[j].satvelxyz[1] = sv_info.velxyz[1];
        PVT_Info[j].satvelxyz[2] = sv_info.velxyz[2];
        PVT_Info[j].satclk_drift = sv_info.clk_drift;
    }else{
        return(retval);
    }

    // build A matrix
    A[j][0] = PVT_Info[j].satvelxyz[0];
    A[j][1] = PVT_Info[j].satvelxyz[1];
    A[j][2] = PVT_Info[j].satvelxyz[2];
    A[j][3] = vector_norm(PVT_Info[j].satposxyz);

}


// **************************************************************
// start of iteration loop

correction = 1e9;
max_iterations = 20;
err_last = 1e9;
for (i=0;i<max_iterations;i++){

    // build the B vector
    for (j=0;j<num_valid_meas;j++){
        tempv[0] = Xc[0] - PVT_Info[j].satposxyz[0];
        tempv[1] = Xc[1] - PVT_Info[j].satposxyz[1];
        tempv[2] = Xc[2] - PVT_Info[j].satposxyz[2];
        Ym = vector_norm(tempv);

        matrix_multiply(1,3,1,PVT_Info[j].satposxyz,PVT_Info[j].satvelxyz,&tempd);

        B[j] = tempd+W[j]*Ym+Xc[3]*(A[j][3]-Ym);
    }


    matrix_transpose(num_valid_meas,4,(double *)A,(double *)Atrans); // At

    matrix_multiply(4,num_valid_meas,4,(double *)Atrans,(double *)A,(double *)AtA);  // AtA
    matrix_multiply(4,num_valid_meas,1,(double *)Atrans,(double *)B,(double *)AtB);  // AtB

    invert4x4(AtA,invAtA);  // inv(AtA)

    // next estimate of receiver values
    matrix_multiply(4,4,1,(double *)invAtA,(double *)AtB,(double *)Xc);   // Xc = inv(A2)*AtB;


//    tempv[0] = Xc[0] - Rx_Pos_Truth[0];
//    tempv[1] = Xc[1] - Rx_Pos_Truth[1];
//    tempv[2] = Xc[2] - Rx_Pos_Truth[2];

    tempv[0] = Xc[0] - Xc_last[0];
    tempv[1] = Xc[1] - Xc_last[1];
    tempv[2] = Xc[2] - Xc_last[2];

    err_last = err;
    err = vector_norm(tempv);

    correction = fabs(err - err_last);

    Xc_last[0] = Xc[0];
    Xc_last[1] = Xc[1];
    Xc_last[2] = Xc[2];

    // if correction is small enough, we're done, exit loop
    if(correction < DOPP_POS_TOL){
        break;
    }

}

// end of iteration loop
// **************************************************************

if (i < max_iterations){

    // success, save result
    misc_info_dopp.recv_pos[0]  = Xc[0];
    misc_info_dopp.recv_pos[1]  = Xc[1];
    misc_info_dopp.recv_pos[2]  = Xc[2];
    misc_info_dopp.clock_drift  = Xc[3];

    system_vars.recv_dopp_pos[0] = Xc[0];
    system_vars.recv_dopp_pos[1] = Xc[1];
    system_vars.recv_dopp_pos[2] = Xc[2];
    system_vars.recv_dopp_pos[3] = Xc[3];

    system_vars.position_status = HAVE_DOPPLER_POS_ESTIMATE;

}else{

    misc_info_dopp.recv_pos[0]  = 99.0;
    misc_info_dopp.recv_pos[1]  = 99.0;
    misc_info_dopp.recv_pos[2]  = 99.0;
    misc_info_dopp.clock_drift  = 0.0;
    // problem, set return flag
    retval = 1;

}

// see how close to truth we got
vector_subtract(Xc,Rx_Pos_Truth,tempv);
misc_info_dopp.final_pos_diff_mag = vector_norm(tempv);
misc_info_dopp.final_iterations = i;

//update_nav_log3();

return(retval);
}

// **************************************************************
//  update_debug_log
// **************************************************************
void update_debug_log3(int num_points)
{
int i;

    if(debug_file3 == NULL){
        debug_file3 = fopen("debug_log3.dat","w");
    }

    if(debug_file3 != NULL){

      if(num_points < 60){

        for(i=0;i<num_points;i++){
        	fprintf(debug_file3, "%15.9f ",misc_info_dopp.debug[i]);
        }
        fprintf(debug_file3, "\n");

      }else{
        fprintf(debug_file3, "XXX\n");
      }

    }
}


// **************************************************************
//  update_nav_log3
// **************************************************************
void update_nav_log3(void)
{

    if(system_vars.nav_log3 == NULL){
        system_vars.nav_log3 = fopen("nav_log3.dat","w");
    }

    if(system_vars.nav_log3 != NULL){

        fprintf(system_vars.nav_log3, "%15.9f ",misc_info_dopp.tg);
        fprintf(system_vars.nav_log3, "%15.9f %15.9f %15.9f ",misc_info_dopp.recv_pos[0],misc_info_dopp.recv_pos[1],misc_info_dopp.recv_pos[2]);
        fprintf(system_vars.nav_log3, "%.5f ",misc_info_dopp.clock_bias);
        fprintf(system_vars.nav_log3, "%.5f ",misc_info_dopp.final_pos_diff_mag);
        fprintf(system_vars.nav_log3, "%.5f ",misc_info_dopp.final_time_diff_mag);
        fprintf(system_vars.nav_log3, "%d ",misc_info_dopp.final_iterations);
        fprintf(system_vars.nav_log3, "\n");

    }
}




