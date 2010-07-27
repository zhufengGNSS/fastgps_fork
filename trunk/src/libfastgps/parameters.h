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

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdio.h>
#include <math.h>
#include "types.h"

#define NAV_OMEGAE_DOT	7.2921151467e-005
#define NAV_C						299792458
#define NAV_GM					3.986005e14
#define NAV_F						-4.442807633e-010
#define NAV_DATA_BIT_DURATION	0.02
#define R2D 57.295779513082320876798154814105           /* radians to degrees */
#define D2R 1.0/R2D                                    /* degrees to radians */

#define NO      0
#define YES     1
#define OFF     0
#define ON      1
#define READY   2
#define OK      1
#define PROBLEM 2
#define SUCCESS 100
#define FAILURE 101

// receiver states
#define LOST            0
#define ACQ             1
#define TRACKING        2
#define NAV             3

#define HAVE_NOTHING    0
#define HAVE_TOW        1
#define HAVE_EPH        2
#define HAVE_EXTERNAL_EPH        1
#define XPVT_INTERVAL    0.2     // seconds
#define NO_LOGGING      0
#define NORMAL_LOGGING  1
#define DEBUG_LOGGING   2
#define GOOGLE_LOGGING  3

#define HAVE_TIME_FILE__ESTIMATE 1
#define HAVE_WGS84_FILE_ESTIMATE 1
#define HAVE_DOPPLER_POS_ESTIMATE 2
#define HAVE_SNAPSHOT_POS_ESTIMATE 3

//#define MAX_CHANNELS 1
#define MAX_CHANNELS 12
#define MAX_SATELLITES 32
#define MAX_SATELLITES_TO_TRACK 8
#define MINIMUM_PVT_SATELLITES 4

#define CHIPS_PER_CODE 1023
#define DATA_BUF_SIZE 2000

#define DOPPLER_RADIUS 5000
#define NUM_COARSE_DOPPLERS 101  // add one for center, keep sides even
#define FINE_DOPPLER_RADIUS (DOPPLER_RADIUS * 2 / NUM_COARSE_DOPPLERS)
#define NUM_FINE_DOPPLERS 251

#define CORRELATOR_BUF_SIZE (2000)
#define CODE_FREQ      1.023e6
#define CARRIER_FREQ   1.57542e9
#define CARRIER_AID_SF (CODE_FREQ / CARRIER_FREQ)

#define ACQ_MS        4
#define COARSE_ACQ_THRESH    25
#define COARSE_ACQ_ADJPEAK  (DOPPLER_RADIUS / NUM_COARSE_DOPPLERS)

// these take a lot of FPU time...
//#define GPS_SIN(x) (sin(x))
//#define GPS_COS(x) (cos(x))

// these macros will use only the sign bit of the trig functions (and thus run faster)
#define GPS_SIN(x) (x > 0 ? 1 : -1)
#define GPS_COS(x) (x > M_PI/2 || x < -M_PI/2 ? -1 : 1)
#define UNWRAP_ANGLE(x) { if (x > M_PI) x -= 2 * M_PI; else if (x < -M_PI) x += 2 * M_PI; }

#define MAKE_LOGS

// These gains calculated "roughly" using the formulas provided in Paul Groves, Artech book
// tracking loop default gains
#define KCO_DEFAULT			0.004 // Kco = 4*BW*CODE_DURATION
#define KCO2_DEFAULT		1.0

// for the first half second, very high BW/Gain to pull the freq in fast
#define KCA2_FLL2_DEFAULT	1.15	// with normalization, BW 50
// for the next half second, drop the loopbandwidth to better determine the freq before jump to PLL 
#define KCA2_FLL1_DEFAULT	0.10	// with normalization, BW 15

#define KCA2_PLL_DEFAULT	0.10	// with normalization, BW 15
#define KCA3_PLL_DEFAULT	0.93	// with normalization, BW 15

// time, in milliseconds, between state transitions
#define TL_FLL_SWITCH_TIME 500
#define TL_PLL_TIME 1000
// must be greater than TL_PLL_TIME
// in the interval between TL_PLL_TIME and TL_PULLIN_TIME the data for detecting the nav bit edge is gathered 
#define TL_PULLIN_TIME 1500		

#define MAX(a,b) (a >= b ? a : b)

#define CH_STATE_UNDEFINED    0
#define CH_STATE_ACQUIRE      1
#define CH_STATE_POSTACQ_SPIN 2
#define CH_STATE_PULLIN       3
#define CH_STATE_TRACKING     4

// shorthand
#define CS_UNDEFINED    CH_STATE_UNDEFINED
#define CS_ACQUIRE      CH_STATE_ACQUIRE
#define CS_POSTACQ_SPIN CH_STATE_POSTACQ_SPIN
#define CS_PULLIN       CH_STATE_PULLIN
#define CS_TRACKING     CH_STATE_TRACKING

#define MS_PER_NAV_BIT 20
#define SUBFRAME_LENGTH 300
#define PAYLOAD_LENGTH 30
#define PREAMBLE_LENGTH 8
#define MAX_PREAMBLE_CANDIDATES 10

#define SF_SEARCHING 0
#define SF_VERIFYING 1
#define SF_WAIT_A_SUBFRAME 2
#define SF_RECEIVING 3

#define TL_MOVE_TO_FINE_TRACKING_COUNT   10
#define TL_MOVE_TO_TRANSITION_COUNT     -1
#define TL_GOOD_DLL_POWER_RATIO          2

#define TWO_NEG_5				0.03125
#define TWO_NEG_29			1.86264514923095703125e-9
#define TWO_NEG_31			4.656612873077392578125e-10
#define TWO_NEG_33			1.16415321826934814453125e-10
#define TWO_NEG_43			1.136868377216160297393798828125e-13
#define TWO_NEG_19			1.9073486328125e-6
#define TWO_NEG_55			2.77555756e-17
#define TWO_POS_4				16

#endif
