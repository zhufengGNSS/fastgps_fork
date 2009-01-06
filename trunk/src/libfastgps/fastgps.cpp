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
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

//////////////////////////////////////////////////////////////////////////////
gps_real_t dopplers[NUM_COARSE_DOPPLERS], fine_dopplers[NUM_FINE_DOPPLERS];
//////////////////////////////////////////////////////////////////////////////

extern FILE *acq_debug;
extern FILE *acq_debug2;
FILE *nav_debug;

// Declare globals
char msg[200];
struct s_system_vars system_vars;

int run_fastgps()
{
  FILE *data_file;
  char *data_buf;
  unsigned i,	sidx;
  unsigned int j;
  unsigned retval;

  fastgps_printf("\r\nWelcome to our software receiver\r\n");

  // initialize misc tracking and nav variables 
  memset(&system_vars,0,sizeof(struct s_system_vars));
  system_vars.Rx_State = ACQ;   // intial state is acquire
  tracking_init();
  nav_init();

  // read the config file
  retval = read_config_file();

  if (retval != OK)
  {

    fastgps_printf("ahhhh couldn't read config file.\nCheck that "
                   "fastgps_config.txt is in executable directory.\n");
    return -1;
  }

  // open the data file and allocate temporary data buffer
  data_file = fopen(system_vars.infilename,"rb");
  data_buf = (char *) malloc(DATA_BUF_SIZE);

  // Set the IGS file name from GPS day and week
  if(system_vars.pvt_aiding_flag)
  {
    fastgps_printf("Using sp3 file for navigation information.\n");
    system_vars.nav_GPS_week = system_vars.datafile_week;
    // STG Test
    // s_SV_Info sv_info;sv_info.prn = 3;sv_info.week = 1463;sv_info.TOW = 385675;
    // retval = GetSVInfo(&sv_info,system_vars.IGSfilename);
  }
  else
  {
      // TODO: determine the week based on the data file
  }

  if (!data_file)
  {
    fastgps_printf("ahhhh couldn't open data file.\nCheck that path to data "
                   "file in fastgps_config.txt is correct.\n");
    return -1;
  }

  // Initialize Doppler bins used during coarse and fine acquisition
  // be warned that the decimal point after the 2's are really necessary!
  // otherwise the compiler generates an unsigned multiply and things get
  // really ugly...
  for (i = 0; i < NUM_COARSE_DOPPLERS; i++)
    dopplers[i] = i * 2.0 * DOPPLER_RADIUS / (NUM_COARSE_DOPPLERS - 1) - 
                  DOPPLER_RADIUS;
  for (i = 0; i < NUM_FINE_DOPPLERS; i++)
    fine_dopplers[i] = i * 2. * FINE_DOPPLER_RADIUS / (NUM_FINE_DOPPLERS - 1) - 
                       FINE_DOPPLER_RADIUS;

  // Initialize tracking channels
  for (i = 0; i < MAX_CHANNELS; i++)
    init_correlator_channel(i);

  // Misc Initialization
  if(system_vars.file_acquisition == OFF)
        init_fft_acq(); // if we are not using an acquisition file

  system_vars.search_channel = 0;
  system_vars.ms_update_flag = OFF;

  double loop_time = DATA_BUF_SIZE/system_vars.sampling_freq;
  double watch_interval = 0.5;     // approximate
  unsigned loops_per_watch = watch_interval/loop_time;
  system_vars.next_pvt_time = 0.0;

  if(system_vars.acq_log_flag == DEBUG_LOGGING)
    fastgps_printf("Generating Acquisition Debug Information.\n");
  if(system_vars.nav_log_flag == DEBUG_LOGGING)
    fastgps_printf("Generating Navigation Debug Information.\n");

  /* ********************* */
  /* Main Loop Starts Here */
  /* ********************* */

  while (!feof(data_file) && (system_vars.process_time < system_vars.run_time))
  {
    // Read Data from file
    size_t bytes_read = fread(data_buf, 1, DATA_BUF_SIZE, data_file);
    unsigned int ch_idx;
    system_vars.process_time = system_vars.loop_count * DATA_BUF_SIZE /
                               system_vars.sampling_freq;
    system_vars.loop_count++;
    // Print message every so often
    if (fmod(system_vars.loop_count,loops_per_watch) == 0)
      fastgps_printf("processing time = %.3f\n", system_vars.process_time);

    // if receiver time has been set, update over this data interval
    if(system_vars.recv_time_valid)
      system_vars.recv_time += (double)bytes_read / system_vars.sampling_freq;

    /* ********************* */
    /*    Acquisition        */
    /* ********************* */

    if(system_vars.Rx_State == ACQ)
    {
      // add more data to the acquisition buffer
      for (sidx = 0; sidx < bytes_read; sidx++)
      {
        acq_buf[acq_buf_write_pos++] = data_buf[sidx];
        if (acq_buf_write_pos >= system_vars.acq_buf_len)
        {
          for (ch_idx = 0; ch_idx < MAX_CHANNELS; ch_idx++)
            c[ch_idx].acq.acq_finish_time = sidx;
          break;  // we have enough data for acquisition
        }
      }
      if(system_vars.file_acquisition == ON)
      {
        // read acquisition data in from file if the acq data buffer is full,
        // perform search, start at the same point, even with file acquisition
        if (acq_buf_write_pos >= system_vars.acq_buf_len)
        {
          retval = read_acquisiton_file();
          if(system_vars.num_channels > 0)
          {
            system_vars.Rx_State = TRACKING;  // we found at least 1 satellite
            fastgps_printf("Satellites acquired using input file.\n"
                           "Starting tracking.\n");
          }
          else
          {
            system_vars.Rx_State = LOST;     // we have found nothing, exit
            fastgps_printf("No satellites acquired using input file.\n");
          }
        }
      }
      else
      {
        // search for satellites specified in the config file
        // if the acq data buffer is full, perform search
        if (acq_buf_write_pos >= system_vars.acq_buf_len)
        {
          // loop through satellites based on config file
          for (j = 0; j < system_vars.num_search_svs; j++)
          {
            // Allocate next satellite to search channel
            c[system_vars.search_channel].prn_num = system_vars.prn_search_list[j];
            fastgps_printf("Searching for PRN %d\n",
                           c[system_vars.search_channel].prn_num);
            // perform search for sat on single channel
            retval = acquire2(data_buf, bytes_read,system_vars.search_channel);
            if(retval == SUCCESS)
            {
              fastgps_printf("Found PRN %d, allocated to channel %d\n", 
                             c[system_vars.search_channel].prn_num,
                             system_vars.search_channel);
              // if satellite found, don't search for it anymore
              system_vars.sats_found |= (1 << c[system_vars.search_channel].prn_num);
              // store channel this sv will be tracked on
              system_vars.prn_to_channel_map[c[system_vars.search_channel].prn_num] = system_vars.search_channel;
              // rotate to next channel
              system_vars.search_channel++;
              init_correlator_channel(system_vars.search_channel);
              system_vars.num_channels++;
            }
            else
              fastgps_printf("Did not find PRN %d\n",
                             c[system_vars.search_channel].prn_num);
            if(system_vars.num_channels >= MAX_CHANNELS)
              break;  // we have sats for all the channels; break out of search
          }  // end of search_svs loop
          // exit acquisition
          if(system_vars.num_channels > 0)
          {
            system_vars.Rx_State = TRACKING;  // we found at least 1 satellite
            fastgps_printf("Starting Tracking.\n");
          }
          else
            system_vars.Rx_State = LOST;     // we have found nothing, exit
          // log results of acquisition search
          update_acq_log();
        }  // end if  acq data buffer is full
      }  // end else file_acquisition
    }     // end of acquisition

    /* ********************* */
    /*    Tracking           */
    /* ********************* */

    if(system_vars.Rx_State >= TRACKING)
    {
      // check if we should use the ephemeris file
      if(system_vars.file_ephemeris == ON)
      {
        retval = read_ephemeris_file();
        system_vars.file_ephemeris++; // only read file once
      }
      system_vars.num_valid_tow = 0;
      system_vars.num_valid_eph = 0;
      for (ch_idx = 0; ch_idx < system_vars.num_channels; ch_idx++)
      {
        c[ch_idx].nav.valid_for_pvt = NO;
        c[ch_idx].nav.pseudorange_valid  = NO;
        if (c[ch_idx].state != CH_STATE_ACQUIRE)
        {
          // if we have allocated an acquired satellite to this channel
          // call correlator and tracking functions
          software_correlator(&c[ch_idx], data_buf + 
                                          c[ch_idx].acq.acq_finish_time,
                              bytes_read - c[ch_idx].acq.acq_finish_time);

          if(c[ch_idx].nav.nav_data_state >= HAVE_TOW)
            system_vars.num_valid_tow++;
          if(c[ch_idx].nav.nav_data_state == HAVE_EPH)
            system_vars.num_valid_eph++;

          // only on the first correlation after acquistion do we adjust for
          // acq_finish_time
          c[ch_idx].acq.acq_finish_time = 0;
        }
      }  // end of channel loop
      // update tracking log
      if ((system_vars.tracking_log_flag == NORMAL_LOGGING && 
           system_vars.ms_update_flag >= 20) ||
          (system_vars.tracking_log_flag == DEBUG_LOGGING && 
           system_vars.ms_update_flag >= 1))
      {
        update_tracking_log();
        system_vars.ms_update_flag = OFF;
      }

      /* ********************* */
      /*    Navigation         */
      /* ********************* */

      // If we have subframe sync and TOW's we can try to perform nav solution
      if (system_vars.num_valid_tow >= 1 && 
          system_vars.process_time > system_vars.next_pvt_time)
      {
        system_vars.next_pvt_time = system_vars.process_time + PVT_INTERVAL;
        system_vars.num_valid_meas = 0;
        // if we have a valid TOW, make a pseudorange measurements
        retval = calc_pseudoranges();
        // loop through the tracking channels, again
        for (ch_idx = 0; ch_idx < system_vars.num_channels; ch_idx++)
        {
          c[ch_idx].nav.valid_for_pvt = NO;
          if(c[ch_idx].nav.pseudorange_valid == YES)
          {
            // calculate satellite positions and clock bias
            if(c[ch_idx].nav.nav_data_state == HAVE_EPH)
            {
              if(!system_vars.recv_time_valid)
              {
                // calculate satellite positions at common receiver time
                retval = ProcessEphemeris(system_vars.nav_GPS_week, 
                    system_vars.nav_GPS_secs,
                    (unsigned int) c[ch_idx].prn_num,
                    &c[ch_idx].nav);
              }
              else
              {
                // calculate satellite positions at transmission time
                retval = ProcessEphemeris(system_vars.nav_GPS_week, 
                    c[ch_idx].nav.tx_time,
                    (unsigned int) c[ch_idx].prn_num,
                    &c[ch_idx].nav);
              }
              if(retval == SUCCESS)
              {
                c[ch_idx].nav.valid_for_pvt = YES;
                system_vars.num_valid_meas++;
              }
            }
            else
            {
              // check for external sat data, i.e. IGS
              if (system_vars.pvt_aiding_flag == HAVE_EXTERNAL_EPH)
              {
                // call IGS routine with week,secs and sv to retrieve sat pos
                // vel and clkbias
                s_SV_Info sv_info;
                sv_info.prn = c[ch_idx].prn_num;
                sv_info.week = system_vars.nav_GPS_week;
                // calc sat positions at common receiver time or tx time
                if (!system_vars.recv_time_valid)
                  sv_info.TOW = system_vars.nav_GPS_secs;
                else
                  sv_info.TOW = c[ch_idx].nav.tx_time;
                retval = GetSVInfo(&sv_info,system_vars.IGSfilename);
                if(retval == 0)
                {
                  c[ch_idx].nav.sat_pos[0] = sv_info.posxyz[0];
                  c[ch_idx].nav.sat_pos[1] = sv_info.posxyz[1];
                  c[ch_idx].nav.sat_pos[2] = sv_info.posxyz[2];
                  c[ch_idx].nav.clock_corr = sv_info.clk_bias;
                  c[ch_idx].nav.sat_vel[0] = sv_info.velxyz[0];
                  c[ch_idx].nav.sat_vel[1] = sv_info.velxyz[1];
                  c[ch_idx].nav.sat_vel[2] = sv_info.velxyz[2];
                  c[ch_idx].nav.clock_drift = sv_info.clk_drift;
                  c[ch_idx].nav.valid_for_pvt = YES;
                  system_vars.num_valid_meas++;
                }
              }
              else
              {
                // no dice, we have to wait for the eph to be decoded in the
                // nav data
              }
            }
            // apply WAAS ionosphere correction
            if (system_vars.Rx_State == NAV && 
                system_vars.waas_flag == YES && 
                c[ch_idx].nav.valid_for_pvt == YES)
              retval = WAAS_corrections(ch_idx);
          }
        } // end if pseudorange OK
      }  // end of channel loop
      // Estimate the receiver position, this assumes everything went perfect
      // above
      if(system_vars.num_valid_meas >= MINIMUM_PVT_SATELLITES)
      {
        retval = PVT_Solution(system_vars.num_channels);
        if(retval == SUCCESS)
          update_nav_log(); // update navigation logs
        if(system_vars.Rx_State != NAV)
        {
          system_vars.Rx_State = NAV;
          fastgps_printf("\nStarting Navigation Calculation.\n\n");
        }
      } // end  Navigation
    }  // end of Rx_State >= TRACKING
  }  // end of while loop

  // clean up and exit
  free(data_buf);
  if(system_vars.file_acquisition == OFF)
    shutdown_fft_acq();

  // close normal log files
  if(system_vars.acq_log != NULL)
    fclose(system_vars.acq_log);
  if(system_vars.acq_log2 != NULL)
    fclose(system_vars.acq_log2);
  if(system_vars.tracking_log != NULL)
    fclose(system_vars.tracking_log);
  if(system_vars.nav_log != NULL)
    fclose(system_vars.nav_log);
  if(system_vars.google_log != NULL)
  {
    /* write the footer to Google Earth file, and close */
    sprintf(msg,"</Folder>\n</kml>");
    fputs(msg,system_vars.google_log);
    fclose(system_vars.google_log);
  }
  // close debug log files
  if(acq_debug != NULL)
    fclose(acq_debug);
  if(acq_debug2 != NULL)
    fclose(acq_debug2);
  if(nav_debug != NULL)
    fclose(nav_debug);
  // close data file
  if(data_file != NULL)
    fclose(data_file);

  fastgps_printf("finished. good night.\n");

  return 0;
}

int read_config_file()
{
  int retval = OK;
  double testd = 0;
  int testi=0;
  char tempc=0;
  int i,num_sats;
  unsigned ch_idx;

  memset(system_vars.infilename,0,200);

  /* Open configuration file */
  system_vars.config_file = fopen("fastgps_config.txt","r");
  if (!system_vars.config_file)
  {
    printf("Couldn't open configuration file 'fastgps_config.txt'.\n");
    exit(1);
  }
  rewind(system_vars.config_file);  // make sure we are at the beginning
  /* Read in info from file */
  if(!system_vars.config_file) 
  {
    retval = 0;
    fastgps_printf("couldn't open config file.\n");
    return retval;
  }
  else
  {
  	while(!feof(system_vars.config_file))  /* until end of file */
    {
      fread(&tempc,1,1,system_vars.config_file);
      if(tempc == '^')
      {
        fread(&tempc,1,1,system_vars.config_file);
        if(tempc == 'A')
        {
          /* Aiding information */
          fscanf(system_vars.config_file," %d",&testi);
          system_vars.pvt_aiding_flag = testi;
          fscanf(system_vars.config_file," %d",&testi);
          system_vars.datafile_week = testi;
          fscanf(system_vars.config_file," %d",&testi);
          system_vars.datafile_day = testi;
          /* IGS Data File */
          fscanf(system_vars.config_file," %s",system_vars.IGSfilename);
        }  // end if 'A'
        if(tempc == 'D') // data file
          fscanf(system_vars.config_file," %s",system_vars.infilename);
        if(tempc == 'E') // file ephemeris flag
        {
          fscanf(system_vars.config_file," %d",&testi);
          system_vars.file_ephemeris = testi;
        }
        if(tempc == 'F')
        {
          fscanf(system_vars.config_file," %lf",&testd);
          system_vars.sampling_freq = testd;
          // other variables dependant on sampling_freq
          system_vars.sample_period = (1.0 / system_vars.sampling_freq);
          system_vars.acq_buf_len = ((int)(ACQ_MS * system_vars.sampling_freq /
                                           1000));
          fscanf(system_vars.config_file," %lf",&testd);
          system_vars.IF = testd;
        }  // end if 'F'
        if(tempc == 'L') // logging flags
        {
          fscanf(system_vars.config_file," %d",&testi);
          system_vars.acq_log_flag = testi;
          fscanf(system_vars.config_file," %d",&testi);
          system_vars.tracking_log_flag = testi;
          fscanf(system_vars.config_file," %d",&testi);
          system_vars.nav_log_flag = testi;
        }  // end if 'L'
        if(tempc == 'S') // satellite info
        {
          fscanf(system_vars.config_file," %d",&testi);
          num_sats = testi;
          if(num_sats == 99) // search all satellites
          {
            for (i = 0; i < MAX_SATELLITES; i++)
              system_vars.prn_search_list[i] = i+1;
            system_vars.num_search_svs = MAX_SATELLITES;
            system_vars.num_sats_to_track = MAX_SATELLITES_TO_TRACK;
          }
          else if(num_sats == 100)
          {
            // use file for acquisition
            system_vars.file_acquisition = ON;
          }
          else if(num_sats < MAX_SATELLITES)
          {
            // read in specific satellties to search for
            for (i = 0; i < num_sats; i++)
            {
              fscanf(system_vars.config_file," %d",&testi);
              system_vars.prn_search_list[i] = testi;
            }
            system_vars.num_search_svs = num_sats;
            system_vars.num_sats_to_track = num_sats;
          }
        }  // end if 'S'
        if(tempc == 'T') // time to process, in seconds
        {
          fscanf(system_vars.config_file," %lf",&testd);
          system_vars.run_time = testd;
        }  // end if 'T'
        if(tempc == 'G')
        {
          int interval1,interval2,interval3;
          double Kco_temp, Kco2_temp, Kca2_fll1_temp, Kca2_fll2_temp;
          double Kca2_pll_temp, Kca3_pll_temp;
          /* Tracking loop intervals and gains */
          fscanf(system_vars.config_file," %d",&testi);
          interval1 = testi;
          fscanf(system_vars.config_file," %d",&testi);
          interval2 = testi;
          fscanf(system_vars.config_file," %d",&testi);
          interval3 = testi;
          fscanf(system_vars.config_file," %lf",&testd);
          Kco_temp = testd;
          fscanf(system_vars.config_file," %lf",&testd);
          Kco2_temp = testd;
          fscanf(system_vars.config_file," %lf",&testd);
          Kca2_fll1_temp = testd;
          fscanf(system_vars.config_file," %lf",&testd);
          Kca2_fll2_temp = testd;
          fscanf(system_vars.config_file," %lf",&testd);
          Kca2_pll_temp = testd;
          fscanf(system_vars.config_file," %lf",&testd);
          Kca3_pll_temp = testd;
          for (ch_idx = 0; ch_idx < MAX_CHANNELS; ch_idx++)
          {
            struct channel *ch = &c[ch_idx];
            ch->track.Kco = Kco_temp;
            ch->track.Kco2 = Kco2_temp;
            ch->track.Kca2_FLL1 = Kca2_fll1_temp;
            ch->track.Kca2_FLL2 = Kca2_fll2_temp;
            ch->track.Kca2_PLL = Kca2_pll_temp;
            ch->track.Kca3_PLL = Kca3_pll_temp;
            ch->track.fll_switch_time = interval1;
            ch->track.pll_switch_time = interval2;
            ch->track.pullin_time = interval3;
          }
        }  // end if 'G'
        if(tempc == 'W') // waas corrections file
        {
          /* WAAS corrections file */
          fscanf(system_vars.config_file," %s",system_vars.WAASfilename);
          system_vars.waas_flag = YES;
        }  // end if 'W'
      }  // if '\'
    }  // end while
  }  // end else
  return retval;
}


void update_acq_log()
{
  if(system_vars.acq_log_flag >= NORMAL_LOGGING)
  {
    // one log used for re-acquisiiton, the other for plotting because
    // Octave/Matlab and Python don't like to load characters into arrays nicer
    // ways to do this are invited.
    if(system_vars.acq_log == NULL)
      system_vars.acq_log = fopen("fastgps_acquisition_log.dat","w");
    if(system_vars.acq_log2 == NULL)
      system_vars.acq_log2 = fopen("fastgps_acquisition_log_plot.dat","w");

    if((system_vars.acq_log != NULL) && (system_vars.acq_log2 != NULL))
    {
      for (int sv = 1; sv <= MAX_SATELLITES; sv++)
      {
        uint64_t tempul64 = (1 << sv);
        if(tempul64 & system_vars.sats_found)
        {
          unsigned tempchan = system_vars.prn_to_channel_map[sv];
          // write acq info for this channel, time and loop reference, start
          // with /A to enable use of this info during future runs with the
          // same data
          fprintf(system_vars.acq_log, "/A %.5f %lu ",
                  system_vars.process_time, system_vars.loop_count);
          fprintf(system_vars.acq_log2, "0 %.5f %lu ",
                  system_vars.process_time, system_vars.loop_count);
          char msg[1024];
          sprintf(msg, "%d %d %d %f %f %f ", sv, tempchan,
                  c[tempchan].prn_num, 0.0, c[tempchan].acq.ratio, 0.0);
          fputs(msg,system_vars.acq_log);
          fputs(msg,system_vars.acq_log2);

          // needed when performing file acquisition
          sprintf(msg, "%12.11f %12.11f %lf %12.5f %12.5f %12.9f %d \n",
                  c[tempchan].true_car_phase_inc,
                  c[tempchan].true_code_inc,
                  c[tempchan].track.carrier_freq,
                  c[tempchan].track.code_freq,
                  c[tempchan].doppler,
                  c[tempchan].code_prompt,
                  c[tempchan].acq.acq_finish_time);
          fputs(msg,system_vars.acq_log);
          fputs(msg,system_vars.acq_log2);
        }
        else
        {
          // write zeros, no sat on this channel, /N indicates satellite was
          // not acquired
          fprintf(system_vars.acq_log, "/N %.5f %lu ",
                  system_vars.process_time, system_vars.loop_count);
          fprintf(system_vars.acq_log2, "0 %.5f %lu ",
                  system_vars.process_time, system_vars.loop_count);
          sprintf(msg, "%d %d %d %f %f %f %f %f %f %f %f %f %d \n",
                  sv, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0);
          fputs(msg,system_vars.acq_log);
          fputs(msg,system_vars.acq_log2);
        }
      }  // end of sv loop
    }
  }
}

void update_tracking_log()
{
  if(system_vars.tracking_log_flag >= NORMAL_LOGGING)
  {
    if(system_vars.tracking_log == NULL)
      system_vars.tracking_log = fopen("fastgps_tracking_log.dat","w");
    if(system_vars.tracking_log != NULL)
    {
      // active channels, time reference
      fprintf(system_vars.tracking_log, "%.5f %d ",
              system_vars.process_time, system_vars.num_channels);
      // loop through channels
      for (unsigned chan = 0; chan < MAX_CHANNELS; chan++)
      {
        fprintf(system_vars.tracking_log, "%d %d %d ",
                c[chan].prn_num, c[chan].ip_save, c[chan].qp_save);
        fprintf(system_vars.tracking_log, "%f ", c[chan].doppler);
        fprintf(system_vars.tracking_log, "%f ",0.0);
      }  // end of channel loop
      fprintf(system_vars.tracking_log, "\n");
    }
  }
}

void update_nav_log()
{
  unsigned chan;
  if(system_vars.nav_log_flag >= NORMAL_LOGGING)
  {
    if(system_vars.nav_log == NULL)
      system_vars.nav_log = fopen("fastgps_navigation_log.dat","w");
    if(system_vars.nav_log != NULL)
    {
      // GPS week and time
      fprintf(system_vars.nav_log, "%.5f %d ",
              system_vars.recv_time, system_vars.num_valid_meas);
      // Estimated PVT info: ECEF, clock bias, etc.
      fprintf(system_vars.nav_log, "%15.5f %15.5f %15.5f %12.9f ",
              system_vars.recv_pos[0], system_vars.recv_pos[1],
              system_vars.recv_pos[2], system_vars.recv_pos[3]);   
      fprintf(system_vars.nav_log, "%12.9f %12.9f %15.5f %12.9f ",
              system_vars.recv_pos_llh[0], system_vars.recv_pos_llh[1],
              system_vars.recv_pos_llh[2], system_vars.gdop);
      fprintf(system_vars.nav_log, "%15.5f %15.5f %15.5f ",
              system_vars.recv_pos_neu[0], system_vars.recv_pos_neu[1],
              system_vars.recv_pos_neu[2]);
      // ECEF velocity, clock drift
      fprintf(system_vars.nav_log, "%15.5f %15.5f %15.5f %12.9f ",
              system_vars.recv_vel[0],system_vars.recv_vel[1],
              system_vars.recv_vel[2],system_vars.recv_vel[3]);
      
      for (chan = 0; chan < MAX_CHANNELS; chan++)
      {
        // 6 entries per channel -> state,PR,satX,satY,satZ,satClkbias
        fprintf(system_vars.nav_log, " %d %15.5f ",
                c[chan].nav.nav_data_state, c[chan].nav.pseudorange);
        fprintf(system_vars.nav_log, "%15.5f %15.5f %15.5f %12.9f ",
                c[chan].nav.sat_pos[0], c[chan].nav.sat_pos[1],
                c[chan].nav.sat_pos[2], c[chan].nav.clock_corr);
        fprintf(system_vars.nav_log, "%15.5f %15.5f %15.5f %12.9f ",
                c[chan].nav.sat_vel[0], c[chan].nav.sat_vel[1],
                c[chan].nav.sat_vel[2],c[chan].nav.clock_drift);
      }
      fprintf(system_vars.nav_log, "\n");
    }
  }

  // Navigation debug info TBD
  if(system_vars.nav_log_flag == DEBUG_LOGGING)
  {
    if(nav_debug == NULL)
      nav_debug = fopen("nav_debug_log.dat","w");
    if(nav_debug != NULL)
      fprintf(nav_debug, "%d %.5f\n", system_vars.num_valid_meas,
              system_vars.nav_GPS_secs); // active channels, time
  }

  if(system_vars.nav_log_flag == GOOGLE_LOGGING)
  {
    if(system_vars.google_log == NULL)
      system_vars.google_log = fopen("fastgps_Google.kml","w");
    if(system_vars.GoogleEarthHeaderWritten == NO)
    {
      sprintf(msg,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
      fputs(msg,system_vars.google_log);
      sprintf(msg,"<kml xmlns=\"http://earth.google.com/kml/2.1\">\n");
      fputs(msg,system_vars.google_log);
      sprintf(msg,"<Folder>\n");
      fputs(msg,system_vars.google_log);
      sprintf(msg,"<Style id=\"fireIcon\">\n");
      fputs(msg,system_vars.google_log);
      sprintf(msg,"<IconStyle>\n<Icon>\n");
      fputs(msg,system_vars.google_log);
      sprintf(msg,"<href>http://maps.google.com/mapfiles/kml/pal3/icon38.png"
                  "</href>\n");
      fputs(msg,system_vars.google_log);
      sprintf(msg,"</Icon>\n</IconStyle>\n</Style>\n");
      fputs(msg,system_vars.google_log);
      system_vars.GoogleEarthHeaderWritten = YES;
    }

    // Add this entry into the GoogleEarth interface file
    sprintf(msg,"<Placemark>\n");
    fputs(msg,system_vars.google_log);
    sprintf(msg,"<name>Rx - GPS sec %ld</name>\n",
            ((unsigned long) system_vars.recv_time));
    fputs(msg,system_vars.google_log);
    sprintf(msg,"<styleUrl>#fireIcon</styleUrl>\n");
    fputs(msg,system_vars.google_log);
    sprintf(msg,"<description>fastgps</description>\n");
    fputs(msg,system_vars.google_log);
    sprintf(msg,"<Point>\n");
    fputs(msg,system_vars.google_log);
    sprintf(msg,"<altitudeMode>absolute</altitudeMode>\n");
    fputs(msg,system_vars.google_log);
    sprintf(msg,"<coordinates>%12.9f,%12.9f,%12.4f</coordinates>\n",
            system_vars.recv_pos_llh[1]*R2D,
            system_vars.recv_pos_llh[0]*R2D,
            system_vars.recv_pos_llh[2]);
    fputs(msg,system_vars.google_log);
    sprintf(msg,"</Point>\n</Placemark>\n");
    fputs(msg,system_vars.google_log);
  }  // if GOOGLE_LOGGING
}

