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
extern "C"
{
  #include "kiss_fft.h"
  #include "kiss_fftr.h"
}

void fine_acquisition(struct channel *ch);

char acq_buf[153000];     // sampling freq 38.192e6 max

unsigned acq_buf_write_pos;

static kiss_fft_cpx *code_fft[32];
static kiss_fft_cfg sample_fft_cfg, inverse_fft_cfg;
static int sample_idx;
static kiss_fft_cpx *sample_buf, *sample_fft, *mult_result, *inverse_fft;
static double *max_shift_energy;
static double  max_doppler, max_shift, max_energy;
extern gps_real_t dopplers[NUM_COARSE_DOPPLERS], 
                  fine_dopplers[NUM_FINE_DOPPLERS];
static double avg_energy;

FILE *acq_debug;
FILE *acq_debug2;

unsigned acquire2(char *samples, unsigned samples_len,unsigned cidx)
{
  unsigned retval = FAILURE;

	// If this channel has already aquired, return
	if (c[cidx].state > CH_STATE_ACQUIRE)
    return(retval);

	if (acq_buf_write_pos >= system_vars.acq_buf_len)
		retval = acq_buffer_full2(cidx);

  return retval;
}

unsigned acq_buffer_full2(unsigned cidx)
{
  double tempd;
  unsigned s; // sample index (saves typing)
  unsigned retval = FAILURE;
  struct channel *ch = &c[cidx];

  for (s = 0; s < system_vars.acq_buf_len; s++)
    max_shift_energy[s] = 0;

  // Doppler loop
  for (ch->acq.doppler_idx = 0; ch->acq.doppler_idx < NUM_COARSE_DOPPLERS; 
       ch->acq.doppler_idx++)
  {
    int nco_sin, nco_cos, i, q;
    ch->car_phase = 0;
    ch->code_prompt = 0;
    avg_energy = 1;
    max_energy = 0;

		ch->car_phase_inc = 2 * M_PI * 
                        (system_vars.IF + dopplers[ch->acq.doppler_idx]) / 
                        system_vars.sampling_freq;

//          sprintf(msg,"Dopp search freq: %d %f \n", ch->acq.doppler_idx,(dopplers[ch->acq.doppler_idx]));
//	        XPRINTF(msg);

		for (s = 0; s < system_vars.acq_buf_len; s++)
		{
			char sample = acq_buf[s];
			nco_sin = GPS_SIN(ch->car_phase);
			nco_cos = GPS_COS(ch->car_phase);
			i =  sample * nco_cos;
			q = -sample * nco_sin;
			ch->car_phase += ch->car_phase_inc;
			UNWRAP_ANGLE(ch->car_phase);
			sample_buf[s].r = i;
			sample_buf[s].i = q;
		}

    kiss_fft(sample_fft_cfg, sample_buf, sample_fft);
		for (s = 0; s < system_vars.acq_buf_len; s++)
		{
			mult_result[s].r = sample_fft[s].r * code_fft[ch->prn_num-1][s].r -
                         sample_fft[s].i * code_fft[ch->prn_num-1][s].i;
			mult_result[s].i = sample_fft[s].r * code_fft[ch->prn_num-1][s].i + 
                         sample_fft[s].i * code_fft[ch->prn_num-1][s].r;
		}
    kiss_fft(inverse_fft_cfg, mult_result, inverse_fft);
		
		// search the result
    for (s = 0; s < system_vars.acq_buf_len; s++)
    {
      double d;
      d = inverse_fft[s].r * inverse_fft[s].r + 
          inverse_fft[s].i * inverse_fft[s].i;
      if (d > max_shift_energy[s % (system_vars.acq_buf_len / ACQ_MS)])
      {
        max_shift_energy[s % (system_vars.acq_buf_len / ACQ_MS)] = d;
        if (d > max_energy)
        {
          max_doppler = dopplers[ch->acq.doppler_idx];
          ch->acq.best_doppler = ch->acq.doppler_idx;
          max_shift = (s % (system_vars.acq_buf_len / ACQ_MS))   / 
                      (double)(system_vars.acq_buf_len / ACQ_MS) * 
                      (double)(CHIPS_PER_CODE);
          max_energy = d;
        }
      }
      avg_energy += d;
    }    

    avg_energy /= (system_vars.acq_buf_len);
    tempd = max_energy / avg_energy;

    if(tempd > COARSE_ACQ_THRESH)
    {
      fastgps_printf("coarse acq best doppler: %f\n", max_doppler);
      // Perform another FFT search at this code delay to narrow the Doppler 
      // frequency to sufficient accuracy to jump right to phase tracking.
      fine_acquisition(ch);

      // Store results for tracking
      ch->true_car_phase_inc = 2 * M_PI * (system_vars.IF + ch->doppler) / 
                               system_vars.sampling_freq;
      ch->true_code_inc = (CODE_FREQ + ch->doppler * CARRIER_AID_SF) / 
                          system_vars.sampling_freq;
      ch->car_phase_inc = 2 * M_PI * (system_vars.IF + 0) / 
                          system_vars.sampling_freq;
      ch->code_inc = (CODE_FREQ + 0 * CARRIER_AID_SF) / 
                     system_vars.sampling_freq;
      ch->track.carrier_freq = system_vars.IF + ch->doppler;
      ch->track.carrier_freq_acq = ch->track.carrier_freq;
      ch->track.code_freq = CODE_FREQ + ch->doppler * CARRIER_AID_SF;
      ch->num_ms = 0;
      ch->tracking_lock = 0;
      ch->ie_accum = ch->qe_accum = 0;
      ch->ip_accum = ch->qp_accum = 0;
      ch->il_accum = ch->ql_accum = 0;
      ch->ie = ch->qe = ch->ip = ch->qp = ch->il = ch->ql = 0;
      ch->state = CH_STATE_POSTACQ_SPIN;
      ch->acq.max_energy = max_energy;
      ch->acq.ratio = tempd;
      retval = SUCCESS;
      fastgps_printf("freq: %f, doppler = %f, detection ratio: %f \n",
                     ch->track.carrier_freq, ch->doppler, tempd);
      break;
    }
  }   // end of Doppler loop

  if(system_vars.acq_log_flag == DEBUG_LOGGING)
  {
    // Acquisition debug information
    if (acq_debug == NULL)
      acq_debug = fopen("acq_debug_log.dat","w");
    if (acq_debug2 == NULL)
      acq_debug2 = fopen("acq_debug2_log.dat","w");

    // STG testing
    complete_delay_search(ch->prn_num,ch->track.carrier_freq);

    // store fft output for Doppler we think signal is at
    for (s = 0; s < system_vars.acq_buf_len; s++)
    {
    	double dd;
    	dd = inverse_fft[s].r * inverse_fft[s].r + 
           inverse_fft[s].i * inverse_fft[s].i;
      fprintf(acq_debug, "%.5f ",dd);
    }
    fprintf(acq_debug, "\n");
  } // end acquisition debug

  return retval;
}

void complete_delay_search(unsigned int prn, double freq)
{
	char nco_sin, nco_cos, *code;
	int code_idx_prompt;
	unsigned s;
  double fa_energy;
  code = CODE_TABLE[prn-1];
  unsigned i;
  double code_step = 0.1;
  unsigned loops = (unsigned) (1023.0/code_step);

  double code_inc = (CODE_FREQ) / system_vars.sampling_freq;
  double car_phase_inc = 2 * M_PI * (freq) / system_vars.sampling_freq;
  double code_prompt = 0;

  for (i = 0; i < loops; i++)
  {
    code_prompt += code_step;
    double car_phase = 0, ip = 0, qp = 0;

    for (s = 0; s < 16368; s++)
    {
      char sample = acq_buf[s];
      nco_sin = sample * GPS_SIN(car_phase);
      nco_cos = sample * GPS_COS(car_phase);
      car_phase += car_phase_inc;
      code_idx_prompt = (int)(code_prompt);
      code_prompt += code_inc;

      ip += code[code_idx_prompt] * nco_sin;
      qp += code[code_idx_prompt] * nco_cos;

      if (car_phase > M_PI)
        car_phase -= 2 * M_PI;
      if (code_prompt > CHIPS_PER_CODE + 1)
        code_prompt -= CHIPS_PER_CODE;
    }
    fa_energy = ip*ip + qp*qp;
    fprintf(acq_debug2, "%.5f ",fa_energy); // store it
  }  // end i loop
  fprintf(acq_debug2, "\n");
}

// **********************************************************************
//  Perform fine Doppler search
// **********************************************************************

void fine_acquisition(struct channel *ch)
{
	char nco_sin, nco_cos, *code;
	int code_idx_prompt;
	unsigned s;
	double fa_max_energy = 0, fa_energy;

  code = CODE_TABLE[ch->prn_num-1];

	for (ch->acq.fine_doppler_idx = 0; 
       ch->acq.fine_doppler_idx < NUM_FINE_DOPPLERS; 
       ch->acq.fine_doppler_idx++)
	{
		ch->doppler = dopplers[ch->acq.best_doppler] + 
                  fine_dopplers[ch->acq.fine_doppler_idx];
		ch->code_inc = (CODE_FREQ + ch->doppler * CARRIER_AID_SF) / 
                   system_vars.sampling_freq;
		ch->car_phase_inc = 2 * M_PI * (system_vars.IF + ch->doppler) / 
                        system_vars.sampling_freq;
		ch->code_prompt = CHIPS_PER_CODE - max_shift + 1;
		ch->ip = ch->qp = 0;

		for (s = 0; s < system_vars.acq_buf_len; s++)
		{
			char sample = acq_buf[s];
			nco_sin = sample * GPS_SIN(ch->car_phase);
			nco_cos = sample * GPS_COS(ch->car_phase);
			ch->car_phase += ch->car_phase_inc;
			code_idx_prompt = (int)(ch->code_prompt);
			ch->code_prompt += ch->code_inc;
			ch->ip += code[code_idx_prompt] * nco_cos;
			ch->qp -= code[code_idx_prompt] * nco_sin;
			if (ch->car_phase > M_PI)
				ch->car_phase -= 2 * M_PI;
			if (ch->code_prompt > CHIPS_PER_CODE + 1)
				ch->code_prompt -= CHIPS_PER_CODE;
		}
		fa_energy = ch->ip*ch->ip + ch->qp*ch->qp;
		if (fa_energy > fa_max_energy)
		{
			fa_max_energy = fa_energy;
			ch->acq.best_fine_doppler = ch->acq.fine_doppler_idx;
		}
	}

	ch->doppler = dopplers[ch->acq.best_doppler] + 
                fine_dopplers[ch->acq.best_fine_doppler];
}

void init_fft_acq()
{
	// sample the C/A code for all satellites
	BYTE sv_num, ch_idx;
	gps_real_t code_time, code_time_inc;
	unsigned s, shift_init;
  kiss_fft_scalar *sampled_code;
  kiss_fftr_cfg forward_fft_cfg;
  //kiss_fft_cpx *fft_buf;

	fastgps_printf("precalculating C/A code FFT tables...");
	sample_idx = 0;
	code_time_inc = 1.0 / system_vars.sampling_freq;
  unsigned acq_fft_len = kiss_fft_next_fast_size(system_vars.acq_buf_len);
	sampled_code = (kiss_fft_scalar *)malloc(sizeof(kiss_fft_cpx) * acq_fft_len);
  memset(sampled_code, 0, acq_fft_len);
  //system_vars.acq_buf_len = 654;
  fastgps_printf("acq buf len = %d\n", acq_fft_len);
  forward_fft_cfg = kiss_fftr_alloc(acq_fft_len, 0, NULL, NULL);
	acq_buf_write_pos = 0;

	for (sv_num = 0; sv_num < MAX_SATELLITES; sv_num++)
	{
		code_fft[sv_num] = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * 
                                              acq_fft_len);
		code_time = 0;
		for (s = 0; s < system_vars.acq_buf_len / ACQ_MS; s++)
		{
			BYTE ms_counter;
			code_time += code_time_inc;
			for (ms_counter = 0; ms_counter < ACQ_MS; ms_counter++)
				sampled_code[s + ms_counter * system_vars.acq_buf_len / ACQ_MS] =
          CODE_TABLE[sv_num][(int)floor(CHIPS_PER_CODE * 1000 * code_time) + 1];
		}
    kiss_fftr(forward_fft_cfg, sampled_code, code_fft[sv_num]);
    // why are we doing this flip? I can't remember
		for (s = 0; s < system_vars.acq_buf_len; s++)
			code_fft[sv_num][s].i *= -1;
	}
	sample_buf = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * acq_fft_len);
	sample_fft = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * acq_fft_len);
	inverse_fft = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * acq_fft_len);
	mult_result = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * acq_fft_len);
  for (unsigned i = 0; i < acq_fft_len; i++)
  {
    sample_buf[i].r = sample_buf[i].i = 0;
    sample_fft[i].r = sample_fft[i].i = 0;
    inverse_fft[i].r = inverse_fft[i].i = 0;
    mult_result[i].r = mult_result[i].i = 0;
  }
  sample_fft_cfg = kiss_fft_alloc(acq_fft_len, 0, NULL, NULL);
  inverse_fft_cfg = kiss_fft_alloc(acq_fft_len, 1, NULL, NULL);
	max_shift_energy = (double *)malloc(sizeof(double) * acq_fft_len);
	for (shift_init = 0; shift_init < acq_fft_len; shift_init++)
		max_shift_energy[shift_init] = 0;
	avg_energy = 0;
	max_energy = 0;
	for (ch_idx = 0; ch_idx < MAX_CHANNELS; ch_idx++)
		c[ch_idx].acq.failure_count = 0;
	fastgps_printf("done\n");
  kiss_fft_free(forward_fft_cfg);
  free(sampled_code);
}

void shutdown_fft_acq()
{
  for (BYTE sv = 0; sv < 32; sv++)
    free(code_fft[sv]);
  kiss_fft_free(sample_fft_cfg);
  kiss_fft_free(inverse_fft_cfg);
  free(sample_buf);
  free(sample_fft);
  free(mult_result);
  free(max_shift_energy);
  free(inverse_fft);
  kiss_fft_cleanup();
}

int read_acquisiton_file()
{
  double testd = 0;
  int testi=0;
  char tempc=0;
  FILE *infile;
  unsigned tempchan = 0;

  /* Open configuration file */
  infile = fopen("fastgps_acquisition_log.dat","r");
  if (!infile)
  {
    printf("woah! couldn't open the acquistion log. please don't select "
           "in the config file.\n");
    exit(1);
  }
  rewind(infile);  // make sure we are at the beginning

  /* Read in info from file */
  if (!infile)
    return system_vars.num_channels;
  else
  {
  	while (!feof(infile))  /* until end of file */
  	{
  		fread(&tempc,1,1,infile);
  		if(tempc == '/')
      {
  	  	fread(&tempc,1,1,infile);
        if(tempc == 'A')
        {
          /* Read SV info entry */
          fscanf(infile," %lf",&testd);   // process time
          fscanf(infile," %d",&testi);
          system_vars.acq_file_start_count = (unsigned) testi;
          fscanf(infile," %d",&testi);  // sv
          fscanf(infile," %d",&testi);  // chan
          tempchan = (unsigned) testi;
          init_correlator_channel(tempchan);
          fscanf(infile," %d",&testi);  // sv again, actual one assigned
          c[tempchan].prn_num = testi;
          fscanf(infile," %lf",&testd);  // max energy
          fscanf(infile," %lf",&testd);   // acq ratio
          fscanf(infile," %lf",&testd);   // Doppler
          // the whole enchilada (for now)
          fscanf(infile," %lf",&testd);
          c[tempchan].true_car_phase_inc = testd;
          fscanf(infile," %lf",&testd);
          c[tempchan].true_code_inc = testd;
          fscanf(infile," %lf",&testd);
          c[tempchan].track.carrier_freq = testd;
          c[tempchan].track.carrier_freq_acq = testd;
          fscanf(infile," %lf",&testd);
          c[tempchan].track.code_freq = testd;
          fscanf(infile," %lf",&testd);
          c[tempchan].doppler = testd;
          fscanf(infile," %lf",&testd);
          c[tempchan].code_prompt = testd;
          fscanf(infile," %lf",&testd);
          c[tempchan].acq.acq_finish_time = testd;
          c[tempchan].state = CH_STATE_POSTACQ_SPIN;
          fastgps_printf("Allocating PRN %d to channel %d\n",
                         c[tempchan].prn_num, tempchan);
          system_vars.num_channels++;
        }  // if 'A'
      }  // if '\'
    }  // end while
  }  // end else
  return system_vars.num_channels;
}

