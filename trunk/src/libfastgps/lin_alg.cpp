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

void invert4x4(double A[4][4], double Ainv[4][4])
{
	double detA = A[0][0]*A[1][1]*A[2][2]*A[3][3]-A[0][0]*A[1][1]*A[2][3]*A[3][2]-A[0][0]*A[2][1]*A[1][2]*A[3][3]+A[0][0]*A[2][1]*A[1][3]*A[3][2]+A[0][0]*A[3][1]*A[1][2]*A[2][3]-A[0][0]*A[3][1]*A[1][3]*A[2][2]-A[1][0]*A[0][1]*A[2][2]*A[3][3]+A[1][0]*A[0][1]*A[2][3]*A[3][2]+A[1][0]*A[2][1]*A[0][2]*A[3][3]-A[1][0]*A[2][1]*A[0][3]*A[3][2]-A[1][0]*A[3][1]*A[0][2]*A[2][3]+A[1][0]*A[3][1]*A[0][3]*A[2][2]+A[2][0]*A[0][1]*A[1][2]*A[3][3]-A[2][0]*A[0][1]*A[1][3]*A[3][2]-A[2][0]*A[1][1]*A[0][2]*A[3][3]+A[2][0]*A[1][1]*A[0][3]*A[3][2]+A[2][0]*A[3][1]*A[0][2]*A[1][3]-A[2][0]*A[3][1]*A[0][3]*A[1][2]-A[3][0]*A[0][1]*A[1][2]*A[2][3]+A[3][0]*A[0][1]*A[1][3]*A[2][2]+A[3][0]*A[1][1]*A[0][2]*A[2][3]-A[3][0]*A[1][1]*A[0][3]*A[2][2]-A[3][0]*A[2][1]*A[0][2]*A[1][3]+A[3][0]*A[2][1]*A[0][3]*A[1][2];

	Ainv[0][0] = (A[1][1]*A[2][2]*A[3][3]-A[1][1]*A[2][3]*A[3][2]-A[2][1]*A[1][2]*A[3][3]+A[2][1]*A[1][3]*A[3][2]+A[3][1]*A[1][2]*A[2][3]-A[3][1]*A[1][3]*A[2][2])/detA;
	Ainv[0][1] = (-A[0][1]*A[2][2]*A[3][3]+A[0][1]*A[2][3]*A[3][2]+A[2][1]*A[0][2]*A[3][3]-A[2][1]*A[0][3]*A[3][2]-A[3][1]*A[0][2]*A[2][3]+A[3][1]*A[0][3]*A[2][2])/detA;
	Ainv[0][2] = (A[0][1]*A[1][2]*A[3][3]-A[0][1]*A[1][3]*A[3][2]-A[1][1]*A[0][2]*A[3][3]+A[1][1]*A[0][3]*A[3][2]+A[3][1]*A[0][2]*A[1][3]-A[3][1]*A[0][3]*A[1][2])/detA;
	Ainv[0][3] = (-A[0][1]*A[1][2]*A[2][3]+A[0][1]*A[1][3]*A[2][2]+A[1][1]*A[0][2]*A[2][3]-A[1][1]*A[0][3]*A[2][2]-A[2][1]*A[0][2]*A[1][3]+A[2][1]*A[0][3]*A[1][2])/detA;

	Ainv[1][0] = (-A[1][0]*A[2][2]*A[3][3]+A[1][0]*A[2][3]*A[3][2]+A[2][0]*A[1][2]*A[3][3]-A[2][0]*A[1][3]*A[3][2]-A[3][0]*A[1][2]*A[2][3]+A[3][0]*A[1][3]*A[2][2])/detA;
	Ainv[1][1] = (A[0][0]*A[2][2]*A[3][3]-A[0][0]*A[2][3]*A[3][2]-A[2][0]*A[0][2]*A[3][3]+A[2][0]*A[0][3]*A[3][2]+A[3][0]*A[0][2]*A[2][3]-A[3][0]*A[0][3]*A[2][2])/detA;
	Ainv[1][2] = (-A[0][0]*A[1][2]*A[3][3]+A[0][0]*A[1][3]*A[3][2]+A[1][0]*A[0][2]*A[3][3]-A[1][0]*A[0][3]*A[3][2]-A[3][0]*A[0][2]*A[1][3]+A[3][0]*A[0][3]*A[1][2])/detA;
	Ainv[1][3] = (A[0][0]*A[1][2]*A[2][3]-A[0][0]*A[1][3]*A[2][2]-A[1][0]*A[0][2]*A[2][3]+A[1][0]*A[0][3]*A[2][2]+A[2][0]*A[0][2]*A[1][3]-A[2][0]*A[0][3]*A[1][2])/detA;

	Ainv[2][0] = (A[1][0]*A[2][1]*A[3][3]-A[1][0]*A[3][1]*A[2][3]-A[2][0]*A[1][1]*A[3][3]+A[2][0]*A[3][1]*A[1][3]+A[3][0]*A[1][1]*A[2][3]-A[3][0]*A[2][1]*A[1][3])/detA;
	Ainv[2][1] = (-A[0][0]*A[2][1]*A[3][3]+A[0][0]*A[3][1]*A[2][3]+A[2][0]*A[0][1]*A[3][3]-A[2][0]*A[0][3]*A[3][1]-A[3][0]*A[0][1]*A[2][3]+A[3][0]*A[0][3]*A[2][1])/detA;
	Ainv[2][2] = (A[0][0]*A[1][1]*A[3][3]-A[0][0]*A[3][1]*A[1][3]-A[1][0]*A[0][1]*A[3][3]+A[1][0]*A[0][3]*A[3][1]+A[3][0]*A[0][1]*A[1][3]-A[3][0]*A[0][3]*A[1][1])/detA;
	Ainv[2][3] = (-A[0][0]*A[1][1]*A[2][3]+A[0][0]*A[2][1]*A[1][3]+A[1][0]*A[0][1]*A[2][3]-A[1][0]*A[0][3]*A[2][1]-A[2][0]*A[0][1]*A[1][3]+A[2][0]*A[0][3]*A[1][1])/detA;

	Ainv[3][0] = (-A[1][0]*A[2][1]*A[3][2]+A[1][0]*A[3][1]*A[2][2]+A[2][0]*A[1][1]*A[3][2]-A[2][0]*A[3][1]*A[1][2]-A[3][0]*A[1][1]*A[2][2]+A[3][0]*A[2][1]*A[1][2])/detA;
	Ainv[3][1] = (A[0][0]*A[2][1]*A[3][2]-A[0][0]*A[3][1]*A[2][2]-A[2][0]*A[0][1]*A[3][2]+A[2][0]*A[0][2]*A[3][1]+A[3][0]*A[0][1]*A[2][2]-A[3][0]*A[0][2]*A[2][1])/detA;
	Ainv[3][2] = (-A[0][0]*A[1][1]*A[3][2]+A[0][0]*A[3][1]*A[1][2]+A[1][0]*A[0][1]*A[3][2]-A[1][0]*A[0][2]*A[3][1]-A[3][0]*A[0][1]*A[1][2]+A[3][0]*A[0][2]*A[1][1])/detA;
	Ainv[3][3] = (A[0][0]*A[1][1]*A[2][2]-A[0][0]*A[2][1]*A[1][2]-A[1][0]*A[0][1]*A[2][2]+A[1][0]*A[0][2]*A[2][1]+A[2][0]*A[0][1]*A[1][2]-A[2][0]*A[0][2]*A[1][1])/detA;

}

double vector_norm(double v[3])
{
    return(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
}

void vector_subtract(double a[3], double b[3],double c[3])
{
    // c = a - b

    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];

}

/********************************************************
* Paramaters:
*           rowsA - number of rows of A
*           columnsA - number of columns of A
*           columnsB - number of columns of B
*           A - pointer to a rowsA x columnsA matrix
*           B - pointer to a coulmnsA x columnsB matrix
*           C - pointer to a rowsA x columnsB matrix
*********************************************************/
void matrix_multiply(unsigned rowsA,unsigned columnsA,unsigned columnsB,double *A,double *B,double *C)
{
unsigned i, j, k;

    for(i = 0;i < rowsA;i++){

        for(j = 0;j < columnsB;j++){

            C[i*columnsB +j] = 0.0;

            for(k = 0;k < columnsA;k++){
                C[i*columnsB +j] +=  A[i*columnsA+k]*B[k*columnsB+j];
            }
        }
    }
}

/********************************************************
* Paramaters:
*           rows - number of rows of A
*           columns - number of columns of A
*           A - pointer to a rows x columns matrix
*           B - pointer to a columns x rows matrix, transpose(A)
*********************************************************/
void matrix_transpose(unsigned rows,unsigned columns,double *A,double *B)
{
unsigned i, j;

    for(i = 0;i < columns;i++){
        for(j = 0;j < rows;j++){
            B[i*rows+j] = A[j*columns+i];
        }
    }
}



