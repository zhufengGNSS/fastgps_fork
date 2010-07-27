/*
* Copyright (c) 2008, Morgan Quigley, Pieter Abbeel and Scott Gleason
* All rights reserved.
*
* Originally written by Morgan Quigley and Pieter Abbeel
* Additional contributions by Scott Gleason
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

#ifndef FUNCS_H
#define FUNCS_H

#include "types.h"

int run_fastgps();

// config and logging
int read_config_file();
void update_acq_log();
void update_tracking_log();
void update_nav_log();

// correlator channel
void init_correlator_channel(uint8_t idx);
void set_state(struct channel *x, int state);
void software_correlator(struct channel *ch, char *samples, unsigned samples_len);

// coarse acquisition
void init_fft_acq();
void shutdown_fft_acq();
void fft_acq_sample(char sample);
void acquire();
unsigned acquire2(char *samples, unsigned samples_len,unsigned chan_num);
void acq_buffer_full();
unsigned acq_buffer_full2(unsigned cidx);
int read_acquisiton_file();
void complete_delay_search(unsigned int prn,double freq);

// tracking
void tracking_init();
void tracking_update(struct channel *x);

// ephemerides
void process_subframe(struct channel *x);
char nav_parity(char *ndat);
void process_nav_bit(struct channel *x, char bit);
void save_ephemeris(struct channel *x);
int read_ephemeris_file();

// navigation
int LeastSquaresPosition(void);
void nav_init();
unsigned int calc_pseudoranges();
int GetSVInfo(struct s_SV_Info *sv, char *infileXname);
unsigned int PVT_Solution(int num_chans);
unsigned int ProcessEphemeris(unsigned int week,double secs, unsigned int sv, nav_info_t *Eph);
unsigned int WAAS_corrections(unsigned int chan);
int ReadWAASFile(double phi_pp_deg,double lambda_pp_deg,double *IGP_points);
double CalculateDelayCorrection(int chan);
void CalculateAzEl(double RxXYZ[3],double SatXYZ[3], double *A,double *E);

// Doppler position
int DopplerPosition(void);
int CalcDopplerPosition(s_nav_packet_dopp *pkt);
void update_debug_log3(int num_points);
void update_nav_log3(void);

// snap shot position
int SingleEpochPosition(void);
int CalcSingleEpochPosition(s_nav_packet *pkt,int signal);
int GetSVInfo2(double *PosVelClk, int prn, int gpsweek, double gpssecs, char *infileXname);
void update_debug_log(int num_points);
void update_nav_log2(void);
double predicted_doppler(int gpsweek, double gpssecs,int prn,double *rxpos,double *rxvel);

// linear algebra
void invert(unsigned actualsize,double *data);
double  vector_dot_product (double *a, double *b);
void invert4x4(double A[4][4], double Ainv[4][4]);
double vector_norm(double v[3]);
void vector_subtract(double a[3], double b[3],double c[3]);
void matrix_multiply(unsigned rowsA,unsigned columnsA,unsigned columnsB, double *A, double *B, double *C);
void matrix_transpose(unsigned rows,unsigned columns,double *A,double *B);

// gnss_utils
void wgsxyz2llh(double *xyz,double *llh);
void wgsllh2xyz(double *llh,double *xyz);
void wgsxyz2neu(double ecef[3],double ref_ecef[3],double NEU[3]);

#endif
