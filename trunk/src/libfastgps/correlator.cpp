/*
* Copyright (c) 2008, Morgan Quigley, Pieter Abbeel and Scott Gleason
* All rights reserved.
*
* Originally written by Morgan Quigley and Pieter Abbeel
* Additional contributions by Scott Gleason
*
* This file is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 2 of the License, or
* (at your option) any later version.
*
* Fastgps is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with fastgps.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "fastgps.h"

//////////////////////////////////////////////////////////////////////////////
struct channel c[MAX_CHANNELS];
//////////////////////////////////////////////////////////////////////////////

void software_correlator(struct channel *ch, char *samples, 
                         unsigned samples_len)
{
	unsigned samples_idx;
	char nco_sin, nco_cos;
	int code_idx_prompt, code_idx_early, code_idx_late;
  char *code;
  char sample;

	code = CODE_TABLE[ch->prn_num-1];
	for (samples_idx = 0; samples_idx < samples_len; samples_idx++)
	{
		sample = samples[samples_idx];
		ch->clock++;

		if (ch->state == CH_STATE_POSTACQ_SPIN)
    {
      // if channel has just acquired, update state and initialize code and phase increments
      // code_prompt has already been set in acquire.cpp
      ch->car_phase_inc = ch->true_car_phase_inc;
      ch->code_inc = ch->true_code_inc;
      ch->state = CH_STATE_PULLIN;
    }

    // advance replica carrier
    nco_sin = sample * GPS_SIN(ch->car_phase);
    nco_cos = sample * GPS_COS(ch->car_phase);
    ch->car_phase += ch->car_phase_inc;

    // advance replica code
    code_idx_prompt = (int)(ch->code_prompt);
    code_idx_early  = (int)(ch->code_prompt + 0.5);
    code_idx_late   = code_idx_early-1;
    ch->code_prompt += ch->code_inc;

    // accumulate at prompt,early and late signals for I and Q
    ch->ip +=  code[code_idx_prompt] * nco_cos;
    ch->ie +=  code[code_idx_early ] * nco_cos;
    ch->il +=  code[code_idx_late  ] * nco_cos;
    ch->qp -=  code[code_idx_prompt] * nco_sin;
    ch->qe -=  code[code_idx_early ] * nco_sin;
    ch->ql -=  code[code_idx_late  ] * nco_sin;

    // adjust for carrier phase roll-overs
    if (ch->car_phase > M_PI)
      ch->car_phase -= 2 * M_PI;

    // if we have processed a complete code phase, update tracking loops
    if (ch->code_prompt > CHIPS_PER_CODE + 1)
    {
      ch->code_prompt -= CHIPS_PER_CODE;

			switch(ch->state)
			{
				case CH_STATE_POSTACQ_SPIN: break;
				case CH_STATE_PULLIN:
				case CH_STATE_TRACKING:  tracking_update(ch); break;
				case CH_STATE_UNDEFINED:
				default: fastgps_printf("ahhhhhhh\n"); break;
			}

      // increment counter for logging of tracking data
      system_vars.ms_update_flag++;
      // save a few things for logging
      ch->ip_save = ch->ip;
      ch->qp_save = ch->qp;

			// reset the code-period accumulators
			ch->ie = ch->qe = ch->ip = ch->qp = ch->il = ch->ql = 0;
		}
	}
}

void init_correlator_channel(uint8_t idx)
{
	struct channel *x = &c[idx];
	x->car_phase = 0;
	x->code_prompt = 0;
	x->code_inc = 0;
	x->sample_idx = 0;
	x->state = CH_STATE_ACQUIRE;
	x->state_ms = 0;
	x->non_coh_disc = 0;
	x->non_coh_disc_integ = 0;
	x->ie_accum = x->qe_accum = x->ip_accum = x->qp_accum = x->il_accum = x->ql_accum = 0;
	x->num_ms = 0;
	x->tracking_lock = 0;
	x->clock = 0;
	x->dll_ctr = 0;
	x->integrated_doppler = 0;
	x->prev_integration_time = 0;
	x->acq.best_doppler = 0;
	x->acq.max_energy = 0;
	x->acq.block_ms = x->acq.block_ctr = 0;
}

void set_state(struct channel *x, int state)
{
	if (state == CH_STATE_ACQUIRE)
	{
		uint8_t ch_idx;

		for (ch_idx = 0; ch_idx < MAX_CHANNELS; ch_idx++)
			if (&c[ch_idx] == x)
				break;

		x->car_phase = 0;
		x->code_prompt = 0;
		x->code_inc = 0;
		x->prn_num = 0;
		x->sample_idx = 0;
		x->state = CH_STATE_ACQUIRE;
		x->acq.best_doppler = 0;
		x->acq.max_energy = 0;
		x->state_ms = 0;
		x->non_coh_disc = 0;
		x->non_coh_disc_integ = 0;
		x->acq.block_ms = x->dll_ctr = x->acq.block_ctr = 0;
		x->ie_accum = x->qe_accum = x->ip_accum = x->qp_accum = x->il_accum = x->ql_accum = 0;
		x->num_ms = 0;
		x->tracking_lock = 0;
		x->clock = 0;
		x->integrated_doppler = 0;
		x->prev_integration_time = 0;
	}
}

