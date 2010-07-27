#ifndef SINGLEEPOCH_H
#define SINGLEEPOCH_H

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


struct s_misc_info
{

    double pos_diff_mag;
    double time_diff_mag;
    double final_pos_diff_mag;
    double final_time_diff_mag;
    unsigned final_iterations;
    double tg_correction;
    double debug[100];
    double pos_correction_mag;

    double recv_pos[3];
    double clock_bias;
    double tg;
    double pdop;

};

struct s_meas
{
	unsigned int    prn;
    double          ms_range;
	double          doppler;

    // secondary values
    double          PR;  // as calculated using the code phase meas
    double          PR2;  // "truth" value

};

struct s_nav_packet
{
	unsigned int    gps_week;
	double          gps_seconds;
	double          time_init;
	double          elapsed_time;
	double          pos_truth[3];
	double          pos_init[3];
	unsigned int    num_meas;

    s_meas          Meas[MAX_CHANNELS];
};

struct s_PVT_Info2
{

   unsigned int         prn;

   double               satposxyz[3];
   double               satvelxyz[3];
   double               doppler;

   double               satclk_bias;
   double               satclk_drift;

   double               TOA_m;
   double               pseudorange;
   double               pseudorange_dot;

};

struct s_globals
{
    int status;
    int nav_count;
    double pos_diff_mag;
    double time_diff_mag;
    double final_pos_diff_mag;
    double final_time_diff_mag;
    int final_iterations;

    double recv_pos[3];
    double recv_pos_llh[3];
    double clock_bias;
    double tg;
    double pdop;
    int nav_position_time_state;
    double pos_correction_mag;
    double tg_correction;
    double correction[5];
    double tg_corr[MAX_CHANNELS];
    double tg_corr_avg;
    double tg_corr_sign[MAX_CHANNELS];
    double tg_corr_sign_best;

    FILE *infile;
    FILE *nav_log;
    FILE *debug_log;
    FILE *google_log;
    int GoogleEarthFileHeader;
    double debug[100];
    double testd;

};

#endif
