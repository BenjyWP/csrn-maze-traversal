#include <math.h>
#include "mex.h"
void maze_neighbor(long maze_size, long cell_num, long direction, long &r, long &c)
// from a given row, col in the maze, returns row and col of
//its neighbor
{
long row, col;
row = (cell_num-1) / maze_size + 1;
col = cell_num - (row-1)*maze_size;
//mexPrintf("row %d col %d\n", row,col);

switch (direction)
{
    case 0: //R
        if (col == maze_size ) {
            c = 1;
            r = row;
        }
        else {
            c = col + 1;
            r = row;
        }
        break;
    case 1: //L
        if (col == 1 ){
            c = maze_size;
            r = row;
        }
        else {
            c = col - 1;
            r = row;
        }
        break;
    case 2: //U
        if (row ==1 ) {
            r = maze_size;
            c = col;
        }
        else {
            r = row - 1;
            c = col;
        }
        break;
    case 3: //D
        if (row ==maze_size ){
            r = 1;
            c = col;
        }
        else {
            r = row + 1;
            c = col;
        }
        break;
}
}

double f(double x)
{
        double z;
        z=(1-exp(-x))/(1+exp(-x));
        return z;
}

/*For Now, Always m=11, N=11, n=5
 REMEMBER -> MATLAB STORES ARRAYS IN COLUMN ORDER !!!
*/
void F_GMLP_cell(/*inputs*/
                double *W, int mW, int nW, double *ww, double *Ws,
                double *x,  double *F_Y, double *F_J_Yhat,
                long n, long m, long N, long cell_num,
                /*outputs*/
                double *F_W, double *F_ww, double *F_Ws, double *F_x )
{
        int i,j;
        //create temp array for F_net
        mxArray *p_F_net;
        double *F_net;
        p_F_net     = mxCreateDoubleMatrix(1, n+m, mxREAL);
        F_net       = mxGetPr(p_F_net);

        for (i=0;i<N;i++) { *(F_x+i)=0; }
        for (i=N;i<N+n;i++) { *(F_x+i)= *(F_Y + i - N); }
        *(F_x + N + n-1) +=  *F_J_Yhat * *Ws;

        //mexPrintf("F_x begin %f\n", *(F_x + N + n-1));

        *F_Ws = *F_Ws + *F_J_Yhat * *(x + N + n-1);

        //mexPrintf("*(x + N + n-1) %f\n", *(x + N + n-1));
        //mexPrintf("F_Ws begin %f\n", *F_Ws);

        *(F_net + N + n - 1)= *(F_x + N + n -1) * (1  -  *(x + N + n - 1) * *(x + N + n - 1))*0.5;

        for (j=0;j<N+n-1;j++)
                *(F_W + mW*j + N + n - 1) =   *(F_net + N + n - 1) * *(x + j);

        for (i=N+n-2;i>=m;i--)
        {
                for (j=i+1;j<N+n;j++)
                {
                        *(F_x+i) += *(W + mW*i + j) * *(F_net + j);
                }
                *(F_net + i)        = *(F_x + i) * (1 - *(x + i) * *(x + i))*0.5;
                *(F_ww  + i - m)    = *(F_x + i) * (1 - *(x + i) * *(x + i))*0.5;
                for (j=0;j<i;j++)
                {
                        *(F_W + mW*j + i) = *(F_net + i) * *(x + j);
                }
        }
        for(i=m-1;i>=0;i--)
                for(j=m;j<N+n;j++)
                       *(F_x + i) += *(W + mW*i + j) * *(F_net + j);

        *(F_ww  +  n - 1)    = *(F_x + N + n -1) * (1 - *(x + N + n -1) * *(x + N + n -1))*0.5;
}

void Maze_Backward_Calc(/*inputs*/
                double *W, int mW, int nW, double *ww, double *Ws, long n, long m,
                long maze_size, long core_iterations,
                double *x,  double *F_Y, double *F_J_Yhat,
                /*outputs*/
                double *F_W, double *F_ww, double *F_Ws, double *F_x, double *store_f_y )
{
    long i, j;
    long c;
    long ro, co;
    long offset, offset2, offset3, offset4, offset5, offset_neigh;
    offset=0; offset2=0; offset3=0; offset4=0; offset5=0;

    for (i=core_iterations; i>0;i--)
    {
       for (c=1;c<=maze_size*maze_size;c++)
       {
           offset  = maze_size*maze_size*(n+m)*(i-1) + (n+m)*(c-1);
           offset2 = maze_size*maze_size*(n+m)*(n+m)*(i-1) + (n+m)*(n+m)*(c-1);
           offset3 = maze_size*maze_size*n*(i-1) + n*(c-1);
           offset4 = (c-1);
           offset5 = (n+m)*(c-1);

            F_GMLP_cell(/*inputs*/ W, mW, nW, ww, Ws, x+offset, F_Y + (c-1)*n, F_J_Yhat + offset4,
                    n, m, m, c, /*outputs*/ F_W+offset2, F_ww+offset3, F_Ws+offset4, F_x +offset5);

       }
       //for (c=1;c<=maze_size*maze_size;c++){mexPrintf("i=%d, F_Ws[%d]= %f\n", i, c, *(F_Ws+c-1));}

       for (j=0; j< n*maze_size*maze_size;j++) *(F_Y+j)=0; // initialization
       for (c=1;c<=maze_size*maze_size;c++)
       {
           offset = (c-1)*(n);
           offset2 = (c-1)*(n+m);
           //neighbor recurrent links
           maze_neighbor(maze_size, c , 3, ro, co); //Down
           offset_neigh =  (m+n) * (maze_size*(ro-1)+ co - 1);
           *(F_Y+offset)+=*(F_x+offset_neigh+2);
           maze_neighbor(maze_size, c , 0, ro, co); //Right
           offset_neigh =  (m+n) * (maze_size*(ro-1)+ co - 1);
           *(F_Y+offset)+=*(F_x+offset_neigh+3);
           maze_neighbor(maze_size, c , 2, ro, co); //Up
           offset_neigh =  (m+n) * (maze_size*(ro-1)+ co - 1);
           *(F_Y+offset)+=*(F_x+offset_neigh+4);
           maze_neighbor(maze_size, c , 1, ro, co); //Left
           offset_neigh =  (m+n) * (maze_size*(ro-1)+ co - 1);
           *(F_Y+offset)+=*(F_x+offset_neigh+5);
           //self recurrent links
           for (j=0; j<n;j++) *(F_Y+offset+j)+=*(F_x+offset2+j+m-n);
       }
       offset2  = maze_size*maze_size*(n)*(i-1);
       //save history
        for (j=0; j<(n)*maze_size*maze_size;j++) *(store_f_y+offset2+j) = *(F_Y+j);
        for (j=0; j<    maze_size*maze_size;j++) *(F_J_Yhat+j)=0;
    }
}

/*[F_W, F_ww, F_Ws, F_x] = BackwardGMLP11_0_5(W, ww, Ws, x, F_Y, F_J_Yhat);*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //inputs
	double *W;
	double *ww;
    double *Ws;
	double *x;
	double *F_Y;
    double *store_f_y;
	double *F_J_Yhat;
    //outputs
    double *F_W;
    double *F_ww;
    double *F_Ws;
    double *F_x;
    //local
    int i, j, ndim;
    int count;
    int nW, mW;
    double n, m, maze_size,core_iterations;

	/*MUST CHECK INPUT FIRST*/
    /*
    for (i=0; i<nrhs; i++)  {
        mexPrintf("\n\n");
        mexPrintf("------------------------------------------------\n");
        mexPrintf("Name: %s%d%c\n", "prhs[",i,']');
        ndim = mxGetNumberOfDimensions(prhs[i]);
        mexPrintf("dimensions %d\n", ndim);
    }
    */
    //get input array pointers
	W           = (double* ) mxGetPr(prhs[0]);
	ww          = (double* ) mxGetPr(prhs[1]);
    Ws          = (double* ) mxGetPr(prhs[2]);
    n           = (long) mxGetScalar(prhs[3]);
    m           = (long) mxGetScalar(prhs[4]);
    maze_size   = (long) mxGetScalar(prhs[5]);
    core_iterations = (long) mxGetScalar(prhs[6]);
    x           = (double* ) mxGetPr(prhs[7]);
    F_J_Yhat    = (double* ) mxGetPr(prhs[8]);

    mW = mxGetM(prhs[0]);
    nW = mxGetN(prhs[0]);

    //create output arrays
	plhs[0] = mxCreateDoubleMatrix(n+m, (n+m)*maze_size*maze_size*core_iterations, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, n*maze_size*maze_size*core_iterations, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, maze_size*maze_size, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, (n+m)*maze_size*maze_size, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, (n)*maze_size*maze_size, mxREAL);
    plhs[5] = mxCreateDoubleMatrix((n)*maze_size*maze_size,core_iterations , mxREAL);

	F_W     = mxGetPr(plhs[0]);
    F_ww    = mxGetPr(plhs[1]);
	F_Ws    = mxGetPr(plhs[2]);
    F_x     = mxGetPr(plhs[3]);
    F_Y     = mxGetPr(plhs[4]);
    store_f_y = mxGetPr(plhs[5]);

    //do the calculations here
    Maze_Backward_Calc(/*inputs*/
                W, mW, nW, ww, Ws, n, m,
                maze_size, core_iterations,
                x,  F_Y, F_J_Yhat,
                /*outputs*/
                F_W, F_ww, F_Ws, F_x, store_f_y );

}
