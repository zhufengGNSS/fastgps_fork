#ifndef DOPP_POS_H
#define DOPP_POS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI                      3.14159265358979
#define TWOPI                   (2.0*PI)
#define RAD2DEG 180/PI
#define DEG2RAD PI/180

#define POS_TOL             1.0    // meters
#define TG_TOL              0.1   // seconds
#define MAX_CHANNELS        12
#define SPEED_OF_LIGHT      299792458.0
#define CODE_FREQ           1.023e6
#define CARRIER_FREQ        1.57542e9
#define L1_WAVELENGTH       (1/CARRIER_FREQ)*SPEED_OF_LIGHT
#define OMEGAE_DOT	        7.2921151467e-005
#define NAV_GM				3.986005e14
#define SECONDS_IN_HOUR     3600
#define SECONDS_IN_DAY      86400
#define R2D 57.295779513082320876798154814105           /* radians to degrees */
#define D2R 1.0/R2D                                    /* degrees to radians */

//#define DOPP_POS_TOL  100.0
#define DOPP_POS_TOL  1.0

struct s_misc_info_dopp
{

    double pos_diff_mag;
    double time_diff_mag;
    double final_pos_diff_mag;
    double final_time_diff_mag;
    int final_iterations;
    double tg_correction;
    double debug[100];
    double pos_correction_mag;

    double recv_pos[3];
    double clock_bias;
    double clock_drift;
    double tg;
    double pdop;

};

struct s_meas_dopp
{
	unsigned int    prn;
	double          doppler;

};

struct s_nav_packet_dopp
{
	unsigned int    gps_week;
	double          gps_seconds;
	double          time_init;
	double          elapsed_time;
	double          pos_truth[3];
	double          pos_init[3];
	unsigned int    num_meas;

    s_meas_dopp     Meas[MAX_CHANNELS];
};

struct s_PVT_Info_dopp
{

   unsigned int         prn;

   double               satposxyz[3];
   double               satvelxyz[3];
   double               doppler;

   double               satclk_bias;
   double               satclk_drift;

};


#endif
