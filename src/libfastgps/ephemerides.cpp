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

static const int preamble[PREAMBLE_LENGTH] = { 1, -1, -1, -1, 1, -1, 1, 1 };

union twobyte_t
{
	short s16;
	unsigned short u16;
} twobyte;

union fourbyte_t
{
	int s32;
	unsigned u32;
} fourbyte;

union onebyte_t
{
	char s8;
	unsigned char u8;
} onebyte;

/* **************************************************************************
* Process Subframe
* ************************************************************************** */

void process_subframe(struct channel *x)
{

  uint8_t *data = x->nav.payload;
  uint8_t subframe_id = (data[5] & 0x1C) >> 2;  // extract subframe id 
  // every subframe contains a TOW, 
  // at this point a psuedorange measurment can be taken
  x->nav.tow = ((data[3] << 9) | (data[4] << 1) | ((data[5] & 0x80) >> 7)) * 6;
  if(x->nav.nav_data_state == HAVE_NOTHING)
    x->nav.nav_data_state = HAVE_TOW;

  fastgps_printf("PRN %d subframe ID = %d  tow = %d\n", 
                 x->prn_num, subframe_id, x->nav.tow);

  // how far are we from the subframe start bit	
  x->nav.subframe_start_clock = x->clock - (unsigned)((NAV_DATA_BIT_DURATION*x->nav.subframe_write_pos)*system_vars.sampling_freq);

  if (!system_vars.recv_time_valid)
  {
		// set the receiver clock
    double time_since_tow = (x->clock - x->nav.subframe_start_clock) / 
                            system_vars.sampling_freq;
    system_vars.recv_time = x->nav.tow + time_since_tow;
    system_vars.recv_time_valid = 1;
		fastgps_printf("Initializing receiver time to %5f \n", 
                   system_vars.recv_time);
  }
  x->nav.subframe_start_valid = 1;
  switch(subframe_id)
  {
    case 1:
      if (!x->nav.subframe_valid[subframe_id])
      {
        x->nav.week_num = (data[6] << 2) + ((data[7] & 0xC0) >> 6);
        onebyte.u8 = data[20];
        x->nav.tgd = onebyte.s8 * TWO_NEG_31;
        x->nav.toc = ((unsigned)(data[22] << 8) + data[23]) * 16;
        onebyte.u8 = data[24];
        x->nav.af2 = ((gps_real_t)onebyte.s8) * TWO_NEG_55;
        twobyte.u16 = (data[25] << 8) + data[26];
        x->nav.af1 = ((gps_real_t)twobyte.s16) * TWO_NEG_43;
        fourbyte.u32 = ((data[27] << 16) + (data[28] << 8) + (data[29] & 0xFC)) << 8;
        fourbyte.s32 >>= 10; // carry the sign bit back down
        x->nav.af0 = ((double)fourbyte.s32) * TWO_NEG_31;
        x->nav.subframe_valid[subframe_id] = 1;
      }
      break;
    case 2:
      if (!x->nav.subframe_valid[subframe_id])
      {
        twobyte.u16 = ((data[7] << 8) + data[8]);
        x->nav.crs = (double)twobyte.s16 * TWO_NEG_5;

        twobyte.u16 = (data[9] << 8) + data[10];
        x->nav.dn = ((gps_real_t)twobyte.s16) * TWO_NEG_43 * M_PI;

        fourbyte.u32 = (data[11] << 24) + (data[12] << 16) + (data[13] << 8) + data[14];
        x->nav.m0 = ((gps_real_t)fourbyte.s32) * TWO_NEG_31 * M_PI;

        twobyte.u16 = (data[15] << 8) + data[16];
        x->nav.cuc = ((gps_real_t)twobyte.s16) * TWO_NEG_29;

        fourbyte.u32 = (data[17] << 24) + (data[18] << 16) + (data[19] << 8) + data[20];
        x->nav.ecc = ((gps_real_t)fourbyte.u32) * TWO_NEG_33;

        twobyte.u16 = (data[21] << 8) + data[22];
        x->nav.cus = ((gps_real_t)twobyte.s16) * TWO_NEG_29;

        fourbyte.u32 = (data[23] << 24) + (data[24] << 16) + (data[25] << 8) + data[26];
        x->nav.sqrta = ((gps_real_t)fourbyte.u32) * TWO_NEG_19;

        x->nav.toe = ((unsigned)(data[27] << 8) + data[28]) * 16;
        x->nav.subframe_valid[subframe_id] = 1;
      }
      break;
    case 3:
      if (!x->nav.subframe_valid[subframe_id])
      {
        twobyte.u16 = (data[6] << 8) + data[7];
        x->nav.cic = twobyte.s16 * TWO_NEG_29;
        fourbyte.u32 = (data[8] << 24) + (data[9] << 16) + (data[10] << 8) + data[11];
        x->nav.omega0 = fourbyte.s32 * TWO_NEG_31 * M_PI;
        twobyte.u16 = (data[12] << 8) + data[13];
        x->nav.cis = twobyte.s16 * TWO_NEG_29;
        fourbyte.u32 = (data[14] << 24) + (data[15] << 16) + (data[16] << 8) + data[17];
        x->nav.inc = fourbyte.s32 * TWO_NEG_31 * M_PI;
        twobyte.u16 = (data[18] << 8) + data[19];
        x->nav.crc = twobyte.s16 * TWO_NEG_5;
        fourbyte.u32 = (data[20] << 24) + (data[21] << 16) + (data[22] << 8) + data[23];
        x->nav.w = fourbyte.s32 * TWO_NEG_31 * M_PI;
        fourbyte.s32 = (data[24] << 24) + (data[25] << 16) + (data[26] << 8);
        fourbyte.s32 >>= 8; // sign-extend it
        x->nav.omegadot = fourbyte.s32 * TWO_NEG_43 * M_PI;
        twobyte.u16 = (data[28] << 8) + data[29];
        twobyte.s16 = twobyte.s16 >> 2;
        x->nav.inc_dot = twobyte.s16 * TWO_NEG_43 * M_PI;
        x->nav.subframe_valid[subframe_id] = 1;
      }
      break;
    case 4:
    case 5:
        // Almanac frames, not used
      break;
  }

  // check if we have all necessary ephemeris subframes, set flag for nav and save them for next time
  if (x->nav.subframe_valid[1] && 
      x->nav.subframe_valid[2] && 
      x->nav.subframe_valid[3])
  {
    if(x->nav.nav_data_state == HAVE_TOW && 
       system_vars.file_ephemeris == OFF)
      save_ephemeris(x);
    x->nav.nav_data_state = HAVE_EPH;
  }

  // STG testing
  //save_ephemeris(x);
}

/* **************************************************************************
*  Process Navigation Data Bit
* ************************************************************************** */

void process_nav_bit(struct channel *x, char bit)
{
  unsigned i;
  char subframe_temp[SUBFRAME_LENGTH];

  // store the bit in the subframe buffer
  x->nav.subframe[x->nav.subframe_write_pos++] = bit;

  if (x->nav.subframe_write_pos >= SUBFRAME_LENGTH)
    x->nav.subframe_write_pos = 0;

  switch(x->nav.subframe_state)
  {
    case SF_SEARCHING:
    {
      // try to correlate with the preamble
      short sf_idx, cor_result = 0;

      // check last 8 bits for preamble pattern
      for (sf_idx = 0; sf_idx < PREAMBLE_LENGTH; sf_idx++)
      {
        cor_result += preamble[sf_idx] * x->nav.subframe[(x->nav.subframe_write_pos - PREAMBLE_LENGTH + sf_idx) % SUBFRAME_LENGTH];
      }
      // if all 8 bits match
      if (abs(cor_result) == PREAMBLE_LENGTH)
      {
        fastgps_printf("found preamble at bit %d\n", x->nav.subframe_write_pos);

        // Store every preamble position found
        if (x->nav.num_preamb_cand < MAX_PREAMBLE_CANDIDATES - 1)
          x->nav.preamb_cand[x->nav.num_preamb_cand++] = 
                                                   x->nav.subframe_write_pos;
      }

      // when we have enough bits for a complete subframe,
      if (x->nav.subframe_write_pos == SUBFRAME_LENGTH - 1)
      {
        uint8_t pc_idx;
        x->nav.subframe_state = SF_VERIFYING;
        fastgps_printf("searching preamble candidates...\n");
        // reset the preamble verification correlation accumulators
        for (pc_idx = 0; pc_idx < MAX_PREAMBLE_CANDIDATES; pc_idx++)
          x->nav.preamb_cor[pc_idx] = 0; 
      }
      break;
    }
    case SF_VERIFYING:
    {
      short sf_idx;
      uint8_t pc_idx;

      // check all the preamble canditates in this 300 bits of data
      for (pc_idx = 0; pc_idx < x->nav.num_preamb_cand; pc_idx++)
      {

        // if the current write position equals that of a previously found preamble,
        //      then check current position for another one.
        //      if there is another one at the current write position
        //          then we have 2 preambles separed by exactly 300 bits
        //              therefore, this is probably a valid subframe start
        //              and the last 300 bits were also a valid subframe

        if (x->nav.subframe_write_pos == x->nav.preamb_cand[pc_idx])
        {
          // test it to see if the preamble appears again at this offset, and that the parity bits check out OK
          for (sf_idx = 0; sf_idx < PREAMBLE_LENGTH; sf_idx++)
            x->nav.preamb_cor[pc_idx] += preamble[sf_idx] * x->nav.subframe[(x->nav.subframe_write_pos - PREAMBLE_LENGTH + sf_idx) % SUBFRAME_LENGTH];

          // if all 8 bits match
          if (abs(x->nav.preamb_cor[pc_idx]) == PREAMBLE_LENGTH)
          {
            char par_buf[32];
            short par_idx;
            char parity;

            // now, build a continuous buffer to check the parity bits
            for (par_idx = 0; par_idx < 32; par_idx++)
              par_buf[par_idx] = x->nav.subframe[(x->nav.subframe_write_pos - PREAMBLE_LENGTH - 2 + par_idx) % SUBFRAME_LENGTH];

            parity = nav_parity(par_buf);

            // if it passes parity check process the subframe, 
			// note: there are cases where this may not be a valid subframe, such as 
			// when the first 8 bits of the HOW are equal to the preamble, but this is OK for now.
            if (parity == 1 || parity == -1)
            {
              unsigned tempi = (x->nav.subframe_write_pos - PREAMBLE_LENGTH) % SUBFRAME_LENGTH;

              fastgps_printf("confirmed preamble at bit %d\n", tempi);

              // process the subframe, then sync for new ones
              // lets re-order the buffer and jump straight to SF_RECEIVING
              unsigned index = 0;
              for (i = tempi; i < SUBFRAME_LENGTH; i++)
                    subframe_temp[index++] = x->nav.subframe[i];
              for (i = 0; i < tempi; i++)
                    subframe_temp[index++] = x->nav.subframe[i];
              for (i = 0; i < SUBFRAME_LENGTH; i++)
                    x->nav.subframe[i] = subframe_temp[i];

              x->nav.parity = parity;
              x->nav.subframe_write_pos = PREAMBLE_LENGTH;
              x->nav.subframe_state = SF_RECEIVING;
              x->nav.first_subframe_flag = 1;

            }  // if parity is OK
          }  // if all 8 bits match
        }  //  if write position == preamble canditate
      }  //  preamble canditate loop

      if (x->nav.subframe_write_pos == SUBFRAME_LENGTH - 1)
      {
        x->nav.subframe_state = SF_SEARCHING; // we didn't find a match... so send it back to the search state
        x->nav.num_preamb_cand = 0;
      }
      break;
    }
    case SF_RECEIVING:
    {
      if (x->nav.subframe_write_pos == 0 || x->nav.first_subframe_flag == 1)
      {
        // this occurs when an entire subframe is valid in the buffer
        short sf_idx, cor_result = 0, word_idx;
        // verify that the preamble is OK
        // check first 8 bits for preamble pattern
        for (sf_idx = 0; sf_idx < PREAMBLE_LENGTH; sf_idx++)
          cor_result += preamble[sf_idx] * x->nav.subframe[sf_idx];

        if (abs(cor_result) == PREAMBLE_LENGTH)
          fastgps_printf("PREAMBLE OK\n");
        else
        {
          // something is wrong, don't process the data
          fastgps_printf("bad preamble.\n");
          x->nav.subframe_state = SF_SEARCHING;
          x->nav.num_preamb_cand = 0;
          break;
        }

        // process the data, skipping the first word
        for (word_idx = 1; word_idx < 10; word_idx++)
        {
          char payload_bits[24];
          // check parity for all words except the first
          for (unsigned pl_idx = 0; pl_idx < 24; pl_idx++)
          {
            payload_bits[pl_idx] = x->nav.subframe[word_idx*30 + pl_idx]; 
            payload_bits[pl_idx] = (payload_bits[pl_idx] > 0 ? 0 : 1);
            char temp_parity = x->nav.subframe[word_idx*30 - 1];  // the last bit in previous word
            if(temp_parity == -1) // invert the bits in this word if needed
              payload_bits[pl_idx] = payload_bits[pl_idx] ^ 1;  // XOR
          }

          // now copy these bits into the payload buffer
          x->nav.payload[word_idx*3] =
            (payload_bits[0] << 7) |
            (payload_bits[1] << 6) |
            (payload_bits[2] << 5) |
            (payload_bits[3] << 4) |
            (payload_bits[4] << 3) |
            (payload_bits[5] << 2) |
            (payload_bits[6] << 1) |
            payload_bits[7];
          x->nav.payload[word_idx*3+1] =
            (payload_bits[ 8] << 7) |
            (payload_bits[ 9] << 6) |
            (payload_bits[10] << 5) |
            (payload_bits[11] << 4) |
            (payload_bits[12] << 3) |
            (payload_bits[13] << 2) |
            (payload_bits[14] << 1) |
            payload_bits[15];
          x->nav.payload[word_idx*3+2] =
            (payload_bits[16] << 7) |
            (payload_bits[17] << 6) |
            (payload_bits[18] << 5) |
            (payload_bits[19] << 4) |
            (payload_bits[20] << 3) |
            (payload_bits[21] << 2) |
            (payload_bits[22] << 1) |
            payload_bits[23];
        }   // end of word loop

        /* **************************** */
        /* Process the subframe data    */
        /* **************************** */
        process_subframe(x);
        x->nav.bit_time = x->nav.tow;
        x->nav.first_subframe_flag = 0;
      }  // if write position == 0
      // update time at this bit, used to take pseudorange measurement
      x->nav.bit_time = x->nav.tow + x->nav.subframe_write_pos * NAV_DATA_BIT_DURATION; // bits come at 50 bps
      break;
    }
  }
}

/* **************************************************************************
*  Check Parity bits
// remember that this function modifies the buffer it is passed!
* ************************************************************************** */
char nav_parity(char *ndat)
{
  uint8_t i;
  char parity[6];
  // check if the data bits must be inverted
  for (i=2; i<24+2; ++i)
    ndat[i] = ndat[1] * ndat[i];

  // calculate parity bits
  parity[0]=ndat[0]*ndat[2]*ndat[3]*ndat[4]*ndat[6]*ndat[7]*
    ndat[11]*ndat[12]*ndat[13]*ndat[14]*ndat[15]*
    ndat[18]*ndat[19]*ndat[21]*ndat[24];

  parity[1]=ndat[1]*ndat[3]*ndat[4]*ndat[5]*ndat[7]*ndat[8]*
    ndat[12]*ndat[13]*ndat[14]*ndat[15]*ndat[16]*
    ndat[19]*ndat[20]*ndat[22]*ndat[25];

  parity[2]=ndat[0]*ndat[2]*ndat[4]*ndat[5]*ndat[6]*ndat[8]*
    ndat[9]*ndat[13]*ndat[14]*ndat[15]*ndat[16]*
    ndat[17]*ndat[20]*ndat[21]*ndat[23];

  parity[3]=ndat[1]*ndat[3]*ndat[5]*ndat[6]*ndat[7]*ndat[9]*
    ndat[10]*ndat[14]*ndat[15]*ndat[16]*ndat[17]*
    ndat[18]*ndat[21]*ndat[22]*ndat[24];

  parity[4]=ndat[1]*ndat[2]*ndat[4]*ndat[6]*ndat[7]*ndat[8]*
    ndat[10]*ndat[11]*ndat[15]*ndat[16]*ndat[17]*
    ndat[18]*ndat[19]*ndat[22]*ndat[23]*ndat[25];

  parity[5]=ndat[0]*ndat[4]*ndat[6]*ndat[7]*ndat[9]*ndat[10]*
    ndat[11]*ndat[12]*ndat[14]*ndat[16]*ndat[20]*
    ndat[23]*ndat[24]*ndat[25];

  // now, check if the computed parity bits equal the received parity bits
  if (parity[0]==ndat[26] && parity[1]==ndat[27] && parity[2]==ndat[28] &&
      parity[3]==ndat[29] && parity[4]==ndat[30] && parity[5]==ndat[31])
    return ndat[1];   // parity is OK. Function output is -1 or 1 depending if the data bits have been inverted (-1) or not (1)
  else
    return 0; // parity check failed. good night.
}

/* **************************************************************************
*  Save Ephemeris
* ************************************************************************** */
void save_ephemeris(struct channel *x)
{
  if(system_vars.nav_log_flag >= NORMAL_LOGGING)
  {
    if(system_vars.eph_log == NULL)
      system_vars.eph_log = fopen("fastgps_ephemeris_log.dat","w");

    if(system_vars.eph_log != NULL)
    {
      // GPS week and time
      fprintf(system_vars.eph_log, "/E %.5f %.5f ",
              system_vars.recv_time, system_vars.process_time);
      fprintf(system_vars.eph_log, " %d %d %d %d ",
              x->prn_num, x->nav.subframe_valid[1], 
              x->nav.subframe_valid[2],x->nav.subframe_valid[3]);
      // subframe 1
      fprintf(system_vars.eph_log, "%d %15.9f %d %15.9f %15.9f %15.9f ",
              x->nav.week_num, x->nav.tgd, x->nav.toc, x->nav.af2,
              x->nav.af1, x->nav.af0);

      // subframe 2
      fprintf(system_vars.eph_log, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %d ",
              x->nav.crs, x->nav.dn, x->nav.m0, x->nav.cuc, x->nav.ecc,
              x->nav.cus, x->nav.sqrta, x->nav.toe);

      // subframe 3
      fprintf(system_vars.eph_log, "%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f ",
              x->nav.cic, x->nav.omega0, x->nav.cis, x->nav.inc, x->nav.crc,
              x->nav.w, x->nav.omegadot, x->nav.inc_dot);

      fprintf(system_vars.eph_log, "\n");
    }
  }
}

/* **************************************************************************/
int read_ephemeris_file()
{
  int num_eph = 0;
  double testd = 0;
  int testi=0;
  char tempc=0;
  FILE *infile;
  unsigned tempchan = 0;
  unsigned found_channel = 0;
  double process_time;
  double eph_time;
  uint8_t temp_prn;

  /* Open configuration file */
  infile = fopen("fastgps_ephemeris_log.dat","r");
  rewind(infile);  // make sure we are at the beginning

  /* Read in info from file */
  if(!infile)
    return(num_eph);
  else
  {
    while (!feof(infile))  /* until end of file */
    {
      VERIFY_IO(fread(&tempc,1,1,infile), 1);
      if(tempc == '/')
      {
        VERIFY_IO(fread(&tempc,1,1,infile), 1);
        if(tempc == 'E')
        {
          /* Read SV info entry */
          VERIFY_IO(fscanf(infile," %lf",&testd), 1); // rcvr_tm of ephemeris 
          eph_time = testd;
          VERIFY_IO(fscanf(infile," %lf",&testd), 1); // process tm since start
          process_time = testd;
          // TODO try setting the rcvr clock
          //			if(!system_vars.recv_time_valid){
          if(0)
          {
            double tdiff = eph_time - process_time;
            system_vars.recv_time = eph_time - tdiff;
            system_vars.recv_time_valid = 1;
            fastgps_printf("Initializing receiver time using eph file to %5f\n",
                           system_vars.recv_time);
          }

          VERIFY_IO(fscanf(infile," %d",&testi), 1);    // prn_num
          temp_prn = (uint8_t) testi;
          VERIFY_IO(fscanf(infile," %d",&testi), 1);    // subframe_valid[0]
          VERIFY_IO(fscanf(infile," %d",&testi), 1);    // subframe_valid[1]
          VERIFY_IO(fscanf(infile," %d",&testi), 1);    // subframe_valid[2]

          // find channel assigned to this prn
          found_channel = 0;
          for (unsigned j = 0; j < system_vars.num_channels; j++)
          {
            if(c[j].prn_num == temp_prn)
            {
              found_channel = 1;
              tempchan = j;
              break;
            }
          }
          if(found_channel)
          {
            // subframe 1
            VERIFY_IO(fscanf(infile," %d",&testi), 1);
            c[tempchan].nav.week_num = testi;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.tgd = testd;
            VERIFY_IO(fscanf(infile," %d",&testi), 1);
            c[tempchan].nav.toc = testi;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.af2 = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.af1 = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.af0 = testd;
            c[tempchan].nav.subframe_valid[1] = 1;

            // subframe 2
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.crs = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.dn = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.m0 = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.cuc = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.ecc = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.cus = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.sqrta = testd;
            VERIFY_IO(fscanf(infile," %d",&testi), 1);
            c[tempchan].nav.toe = testi;
            c[tempchan].nav.subframe_valid[2] = 1;

            // subframe 3
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.cic = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.omega0 = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.cis = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.inc = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.crc = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.w = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.omegadot = testd;
            VERIFY_IO(fscanf(infile," %lf",&testd), 1);
            c[tempchan].nav.inc_dot = testd;
            c[tempchan].nav.subframe_valid[3] = 1;

            fastgps_printf("Have ephemeris for PRN %d on channel %d\n",
                           temp_prn, tempchan);

            num_eph++;
          }  // if found_channel
        }  // if 'E'
      }  // if '\'
    }  // end while
  }  // end else
  return num_eph;
}

