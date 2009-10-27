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

#include "fastgps.h"
#include "types.h"

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void tracking_update(struct channel *x)
{
  gps_real_t E, L, P;
  double Kca2_fll;     // FLL constant
  double diff_codephase_error;
  double diff_carrier_error;

    // *********************************
    // Early, Late and Prompt magnitudes
    // *********************************

	E = sqrt((double)x->ie * x->ie + x->qe * x->qe);
	L = sqrt((double)x->il * x->il + x->ql * x->ql);
	P = sqrt((double)x->ip * x->ip + x->qp * x->qp);

	x->track.P_filt = (9*x->track.P_filt + P)/10;  // simple filter

    // ************************
    // Bit Sync
    // ************************

	if (x->state_ms < x->track.pullin_time && x->state_ms > x->track.pll_switch_time){

        // before we start decoding nav data,
        // find the start of the nav bit during the pullin interval,
        // by watching for the bit flip on the I prompt correlator
		x->track.nav_deriv[x->state_ms % MS_PER_NAV_BIT] += abs(x->ip - x->track.ip_last);

		x->track.code_lock_threshold = x->track.P_filt/2;  // if signal magnitude drops by half we've lost signal

    }else if (x->state_ms == x->track.pullin_time){

        // mark the nav bit start ms from among the 20 possibilities (i.e. 20ms per L1 C/A code navigation data bit)
        // this is important for tracking the signal transmission time and making pseudorange measurements
		uint8_t i, max_nav_idx = 0;
		gps_real_t max_nav = 0;
		for (i = 0; i < MS_PER_NAV_BIT; i++)
			if (x->track.nav_deriv[i] > max_nav)
			{
				max_nav_idx = i;
				max_nav = x->track.nav_deriv[i];
			}
		x->track.nav_bit_start = max_nav_idx;
	}

    // accumulate I prompt over 20ms ...
	if (x->state_ms >= x->track.pullin_time)
	{
		x->track.nav_ctr++;
		x->track.nav_sum += x->ip;
	}

    // ************************
    // Process Nav Bits
    // ************************

	// then collect the nav bit every 20ms
	if ((x->state_ms - x->track.nav_bit_start) % MS_PER_NAV_BIT == 0)
	{
		if (x->track.nav_ctr == MS_PER_NAV_BIT)
		{
			// nav bit is determined using accumulated I prompt values
			int nav_bit = (x->track.nav_sum > 0 ? 1000 : -1000);

            // save nav bit
            x->track.current_nav_bit = (nav_bit > 0 ? 1 : -1);

            /* process each nav bit as its received */
            process_nav_bit(x, x->track.current_nav_bit);

		}
		// reset bit counter and sum
		x->track.nav_ctr = 0;
		x->track.nav_sum = 0;
	}

    // ********************************
	// Tracking loops
    // ********************************

	// calculate code error discriminator, and normalize it
	x->track.code_error = (E - L) / (P + P);  // disc/normalization

	if (x->track.carr_mode){

		/* *********************** */
		// PLL
		/* *********************** */

        if(x->ip == 0){
       		x->track.carr_error = 0.0;
        }else{
       		// calculate phase error discriminator, not normalized by (2*PI)
			// (2*PI) absorbed into gains Kca2, Kca3 
       		x->track.carr_error = atan( (double)x->qp / (double) x->ip);  // radians
        }


		// print message when we have reached the PLL
		if (x->track.carr_mode == 1)
    {
			x->track.old_carr_error = x->track.carr_error;
			x->track.carr_mode = 2;
			x->track.P_filt = P;
			fastgps_printf("reached phase lock loop\n");
		}

		// calculate change in carrier error
		diff_carrier_error = (x->track.carr_error - x->track.old_carr_error);
		// update carrier nco
    x->track.carr_nco = x->track.old_carr_nco + (x->track.Kca2_PLL)*x->track.carr_error + (x->track.Kca3_PLL)*diff_carrier_error;
	}
  else
  {
		/* *********************** */
		// FLL
		/* *********************** */
		x->track.carr_error = 0.0;
		if (x->state_ms >= 1)
    {
      // frequency descriminator, 4 quadrant arctan, atan(cross,dot)
      double num = x->qp*x->track.ip_last - x->ip*x->track.qp_last;  
      double den = x->ip*x->track.ip_last + x->qp*x->track.qp_last;  

      if(den == 0)
        x->track.carr_error = 0.0;
      else
        // calculate freq error discriminator, not normalized by (2*PI)
        x->track.carr_error = atan2(num,den);
    }

		// FLL gain switch
		if (x->state_ms >= x->track.fll_switch_time)
			// after a short interval, drop the loopbandwidth to better determine the freq before jump to PLL 
			Kca2_fll = x->track.Kca2_FLL1;
		else
			// for the first half second, very high BW/Gain to pull the freq in fast
			Kca2_fll = x->track.Kca2_FLL2;     
		
		// update carrier nco using only first order error
		x->track.carr_nco = x->track.old_carr_nco + (Kca2_fll)*x->track.carr_error;

		// timed switch to PLL
		// after pull in time (1.5 sec), switch to PLL
		if (x->state_ms >= x->track.pll_switch_time){
			x->track.carr_mode = 1;
		}


	}   // end else
	

	/* *********************** */
	// DLL
	/* *********************** */
	
	// calculate change in code error
	diff_codephase_error = (x->track.code_error - x->track.old_code_error);
	// update code nco
	x->track.code_nco = x->track.old_code_nco + (x->track.Kco)*x->track.code_error + (x->track.Kco2)*diff_codephase_error;
	
	// ip and qp for next time
	x->track.ip_last = x->ip;
	x->track.qp_last = x->qp;

	// loss of code lock detector
	if(x->state_ms >= x->track.pullin_time && x->track.P_filt < x->track.code_lock_threshold)
  {
		fastgps_printf("ahhh lost signal code lock, prn %d, mag %.1f\n",
                   x->prn_num, P);
		set_state(x, CH_STATE_ACQUIRE);
		//return;
	}

	/* *********************** */
	// update correlators
	/* *********************** */

	if (x->state_ms >= 1)
  {
    x->track.old_carr_nco = x->track.carr_nco;
    x->track.old_carr_error = x->track.carr_error;
    x->track.carrier_freq = x->track.carrier_freq_acq + x->track.carr_nco;
    x->doppler =  x->track.carrier_freq_acq + x->track.carr_nco - system_vars.IF;

    x->track.old_code_nco = x->track.code_nco;
    x->track.old_code_error = x->track.code_error;
    x->track.code_freq = CODE_FREQ + x->track.code_nco;
    x->car_phase_inc = 2 * M_PI * (x->track.carrier_freq) / system_vars.sampling_freq;
    x->code_inc = (x->track.code_freq) / system_vars.sampling_freq;
	}

  // keep track of the number of codes processed
  x->state_ms++;
}

void tracking_init()
{
	unsigned ch_idx;
	for (ch_idx = 0; ch_idx < MAX_CHANNELS; ch_idx++)
	{
		struct channel *ch = &c[ch_idx];
		ch->track.nav_bit_start = 0;
		ch->track.nav_sum = 0;
		ch->track.need_qel_accums = 1;

		ch->track.code_error = ch->track.old_code_error = 0;
		ch->track.carr_error = ch->track.old_carr_error = 0;
		ch->track.code_nco = ch->track.old_code_nco = 0;
		ch->track.carr_nco = ch->track.old_carr_nco = 0;	
		ch->track.ip_last = ch->track.qp_last = 0;
		ch->track.carr_mode = 0;

		// set default tracking gains, these can be changed by config file command
		ch->track.Kco = KCO_DEFAULT;
		ch->track.Kco2 = KCO2_DEFAULT;
		ch->track.Kca2_FLL1 = KCA2_FLL1_DEFAULT;
		ch->track.Kca2_FLL2 = KCA2_FLL2_DEFAULT;
		ch->track.Kca2_PLL = KCA2_PLL_DEFAULT;
		ch->track.Kca3_PLL = KCA3_PLL_DEFAULT;
		
		// default loop switching times
		ch->track.fll_switch_time = TL_FLL_SWITCH_TIME;
		ch->track.pll_switch_time = TL_PLL_TIME;
		ch->track.pullin_time = TL_PULLIN_TIME;

	}

}


