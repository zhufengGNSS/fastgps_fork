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

#ifndef STRUCTS_H
#define STRUCTS_H

#include "parameters.h"

struct tracking_coeffs_t
{
	gps_real_t pll_kp, pll_ki;
	gps_real_t dll_kp, dll_ki;
	uint8_t pll_int_time, dll_int_time, disc_type;
};

struct tracking_info_t
{
	gps_real_t carrier_freq, carrier_freq_acq, code_freq;
	unsigned nav_bit_start, nav_ctr;
	gps_real_t nav_sum;
	gps_real_t nav_deriv[MS_PER_NAV_BIT];
	char need_qel_accums;
    int current_nav_bit;
		
	gps_real_t code_error,old_code_error,carr_error,old_carr_error;
	gps_real_t code_nco,old_code_nco,carr_nco,old_carr_nco;	
	gps_real_t Kco,Kco2,Kca2_FLL1,Kca2_FLL2,Kca2_PLL,Kca3_PLL;
	int ip_last,qp_last;
	int carr_mode;
	unsigned fll_switch_time, pll_switch_time, pullin_time;
	gps_real_t code_lock_threshold,P_filt;

};

struct nav_info_t
{
	char subframe[SUBFRAME_LENGTH];
	uint8_t payload[PAYLOAD_LENGTH];
	unsigned subframe_write_pos, preamb_cand[MAX_PREAMBLE_CANDIDATES];
	uint8_t subframe_state, num_preamb_cand, first_subframe_flag;
	char parity;
	char preamb_cor[MAX_PREAMBLE_CANDIDATES];
	unsigned tow, week_num;
	gps_real_t crs, crc, cuc, cus, cic, cis;
	gps_real_t dn, m0, ecc, sqrta, omega0, omegadot, w, inc, inc_dot;
	gps_real_t af0, af1, af2;
	unsigned toe, toc;
	double tgd;
	uint8_t subframe_valid[6], subframe_start_valid;
	unsigned subframe_start_clock;
	double sat_pos[3], clock_corr, sat_vel[3],clock_drift, pseudorange, old_pseudorange;
	unsigned pseudorange_valid, pseudorange_count;
	unsigned valid_for_pvt;
	double recv_pos[4];
	double bit_time,tx_time;
  uint8_t nav_data_state;
	double doppler_meas;
};

struct acq_info_t
{
	uint8_t state;
	uint8_t doppler_idx, best_doppler;
	uint8_t fine_doppler_idx, best_fine_doppler;
	unsigned short block_ms, block_ctr;
	gps_real_t energy, max_energy;
	gps_real_t ratio;
	unsigned failure_count;
	unsigned acq_finish_time;
};

struct channel
{
	uint8_t prn_num;
	unsigned sample_idx;
	int ie, qe, ip, qp, il, ql, ire, qre, irl, qrl;
	int ip_save, qp_save;
	int ie_accum, qe_accum, ip_accum, qp_accum, il_accum, ql_accum, ire_accum, qre_accum, irl_accum, qrl_accum;
	gps_real_t car_phase, car_phase_inc, true_car_phase_inc;
	gps_real_t code_prompt, code_inc, true_code_inc; // code phase E,L are generated from this on-the-fly
	gps_real_t non_coh_disc, non_coh_disc_integ, doppler;
	uint8_t state;
	unsigned num_ms, state_ms, dll_ctr;
	unsigned clock, prev_integration_time;
	int tracking_lock;
	gps_real_t integrated_doppler; // aka integrated carrier phase
	struct acq_info_t acq;
	struct tracking_info_t track;
	struct nav_info_t nav;
        gps_real_t tx_time;
};


struct s_SV_Info
{

   unsigned int         prn;
   double               posxyz[3];
   double               velxyz[3];

   double               clk_bias;
   double               clk_drift;

   unsigned int         week;
   double               TOW;

};

struct s_PVT_Info
{

   double               satposxyz[3];
   double               satvelxyz[3];
   double               doppler;

   double               satclk_bias;
   double               satclk_drift;

   double               pseudorange;
   double               pseudorange_dot;

};


struct s_system_vars
{

        int     Rx_State;

        unsigned int datafile_week;
        unsigned int datafile_day;

        gps_real_t  sampling_freq;
        gps_real_t  IF;
        gps_real_t  sample_period;
        unsigned int acq_buf_len;

        // Acquisition info
        unsigned int prn_search_list[MAX_SATELLITES];
        unsigned int num_search_svs;
        unsigned int num_channels;
        uint64_t sats_found;
        unsigned int prn_to_channel_map[MAX_SATELLITES];
        unsigned file_acquisition;
        unsigned acq_file_start_count;

        double run_time;
        unsigned long current_ms;
        unsigned long loop_count;
        double process_time;

        int acq_log_flag;
        int tracking_log_flag;
        int nav_log_flag;

        // tracking channels
//        struct channel chan[MAX_CHANNELS];
        unsigned int search_channel;
        unsigned int num_sats_to_track;
        unsigned int ms_update_flag;

        // navigation information
        unsigned file_ephemeris;
       	double initial_pos_guess[3];  // very rough guess of receiver location
       	double recv_pos[4];  // extra value for clock bias
       	double gdop;
       	double recv_pos_refxyz[3];
       	double recv_pos_llh[3];
       	double recv_pos_neu[3];
        double recv_vel[4];     // extra value for clock drift
        unsigned last_nav_time;
        unsigned num_valid_meas;
        unsigned pvt_aiding_flag;
        unsigned waas_flag;
        unsigned num_valid_tow;
        unsigned num_valid_eph;
        double next_pvt_time;
        unsigned nav_GPS_week;
        double nav_GPS_secs;
        double recv_time;
        unsigned recv_time_valid;

        // doppler and snap shot position 
        int snap_shot_flag;
        int doppler_pos_flag;
        int estimate_gpsweek;
        double estimate_gpssecs;
        double estimate_wgs84_pos[3];
        int time_status;
        int position_status;
       	double recv_dopp_pos[4];  // coarse estimate from doppler, extra value for clock drift
       	double recv_snapshot_pos[4];  // snap shot position estimate, extra value for clock bias
       	double recv_snapshot_time;  // snap shot time estimate

        // files
        unsigned GoogleEarthHeaderWritten;
        char infilename[200];
        char WAASfilename[200];
        char IGSfilename[200];
        FILE *config_file;
        FILE *data_file;
        FILE *waas_igp_file;

        FILE *acq_log;
        FILE *acq_log2;
        FILE *tracking_log;
        FILE *tracking_log_test;
        FILE *IQ_log;
        FILE *nav_log;
        FILE *nav_log2;
        FILE *nav_log3;
        FILE *eph_log;
        FILE *google_log;

        double PVT_INTERVAL;
        double LOOP_TIME;
        int testi;
};

#endif
