#include <math.h>
#include "mex.h"


/*
    isParetoSetMember returns the logical Pareto membership of a set of points.

    synopsis:   membership = isParetoSetMember(objectiveMatrix)

    powered by Gianluca Dorini
    
    g.dorini@ex.ac.uk
    
    for compiling type 

    mex -DranSHR3 isparetosetMember.c
 *
 *% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.


    
*/


void isparetosetMember(double * answer, double * M, unsigned int row, unsigned int col);

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
	char * answer;
	double * M;
	unsigned int row, col;
	const int  *dims;
    
	if(nrhs == 0 || nlhs > 1)
	{
	    printf("\nsynopsis:   membership = isParetoSetMember(objectiveMatrix)");
	    plhs[0]    = mxCreateDoubleMatrix(0 , 0 ,  mxREAL);
	    return;
	}
	
	M = mxGetPr(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	row = dims[0];
	col = dims[1];
	
	
	
	/* ----- output ----- */

	plhs[0]    = mxCreateDoubleMatrix(row , 1 ,  mxREAL);
	answer          = mxGetPr(plhs[0]);
	
	
	/* main call */
	isparetosetMember(answer,  M, row, col);
	

}


void isparetosetMember(double * answer, double * M, unsigned int row, unsigned int col)
{
    unsigned int t,i,j;
    double *P, sum;
    bool dominated;
    
    
    P = (double *)mxMalloc(col*sizeof(double));
    for(t = 0; t<row; t++) answer[t] = 1;
    
    for(t = 0; t<row; t++)
    {
        for(j = 0; j<col;j++) P[j] = M[t + j*row]; 
        i = 0;
        dominated = false;
        while(i<row && !dominated)
        {
            if(answer[i])
            {
                j = 0;
                sum = 0;
                while(j < col && M[i + j*row] <= P[j])
                {
                    sum += M[i + j*row] - P[j];
                    j++;
                }
                if(j == col && sum != 0) 
                {
                    answer[t] = false;
                    dominated = true;
                }
            }
            i++;
        }
        
    
    
    }
 
    mxFree(P);
}