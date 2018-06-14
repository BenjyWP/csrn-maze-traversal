#include <stdlib.h>
#include <math.h>
#include "mex.h"

void maze_neighbor(long maze_size, long cell_num, long direction, long &r, long &c)
// from a given row, col in the maze, returns row and col of
//its neighbor
{
long row, col;
row = (cell_num-1) / maze_size + 1;
col = cell_num -  (row-1)*maze_size;
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
OUTPUTS ARE
	y[0] .. y[4]
	v[0] .. v[4]
*/
void GMLP_cell( double *W, long mW, long nW, double *ww,double *xx,
                double *y, long n, long m, long N, double *v, long cell_num)
{
        long i,j;
        long offset, offset2;
        double net;

        offset= (cell_num-1)*(n+m);
        offset2= (cell_num-1)*(n);
        //mexPrintf("current start %d , %d\n", offset, offset2);

        //copy input into local array
        for (i=0; i<m;i++) *(y+offset+i) = *(xx+offset+i);
        for (i=m; i<N+n;i++) *(y+offset+i) = 0;

        for (i=m;i<N+n;i++)
        {
                net=0;
                for (j=0;j<i;j++)
                {
                        net += *(W + mW*j+i) * *(y+offset+j);
                }
                *(v + offset2 + i - m) =   net + *(ww + i-m);
                *(y + offset + i)     = f(net + *(ww + i-m));
        }
}

/*x is the total input to the maze, 16*49 by 1, and y is the total output*/
void Maze_Forward_Calc(/*inputs*/double *W, long mW, long nW, double *ww, double *x
                                ,long n, long m,  long maze_size, long core_iterations,
                        /*outputs*/ double *y, double *v, double *store_y,
                        /*test args*/long test_cell_num, long test_source, long test_target, double* test_dw)
{
    long i, j;
    long c,offset, offset_neigh;
    long ro, co;
    long store_offset;

    for (i=1; i<=core_iterations;i++)
    {
       for (c=1;c<=maze_size*maze_size;c++)
       {
           //test changes one of the weights, if test_cell_num==-1, in all
           //"columns", otherwise in only one column test_cell_num
           if ((test_cell_num == c) || (test_cell_num==-1)){
               //mexPrintf("Adjusted Weight is %f\n", *(W + mW*(test_target-1)+test_source-1));
               *(W + mW*(test_target-1)+test_source-1)+=*test_dw;
               GMLP_cell(W, mW, nW, ww, x, y, n,m, m,  v, c);
               *(W + mW*(test_target-1)+test_source-1)-=*test_dw;
           }
           else {
               GMLP_cell(W, mW, nW, ww, x, y, n,m, m,  v, c);
           }
       }
       //store values
       store_offset = (i-1)*(m+n)*maze_size*maze_size;
       //mexPrintf("store offset %d\n", store_offset);
       for (j=0; j<(m+n)*maze_size*maze_size;j++) *(store_y+store_offset+j) = *(y+j);

       //recurrent links now
       for (c=1;c<=maze_size*maze_size;c++)
       {
           offset = (c-1)*(n+m);
           //mexPrintf("recc offset %d -> %d\n", c, offset);
           //self recurrent
           for (j=m; j<m+n;j++) *(x+offset+j-n) = *(y+offset+j);
           //neighbors

           maze_neighbor(maze_size, c , 3, ro, co); //Down
           offset_neigh =  (m+n) * (maze_size*(ro-1)+ co - 1);
           *(x+offset_neigh+2) = *(y+offset+m);
           maze_neighbor(maze_size, c , 0, ro, co); //Right
           offset_neigh =  (m+n) * (maze_size*(ro-1)+ co - 1);
           *(x+offset_neigh+3) = *(y+offset+m);
           maze_neighbor(maze_size, c , 2, ro, co); //Up
           offset_neigh =  (m+n) * (maze_size*(ro-1)+ co - 1);
           *(x+offset_neigh+4) = *(y+offset+m);
           maze_neighbor(maze_size, c , 1, ro, co); //Left
           offset_neigh =  (m+n) * (maze_size*(ro-1)+ co - 1);
           *(x+offset_neigh+5) = *(y+offset+m);

       }
    }
}

/*[y, v]=ForwardGMLP11_0_5(W, ww, x)*/
// inputs prhs[0] -> W
// input prhs[1] -> ww
// input prhs[2] -> x
// input prhs[3] -> n
// input prhs[4] -> m
// input prhs[5] -> maze_size
// input prhs[6] ->core_iterations
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *W;
	double *ww;
	double *x;
	double *y;
	double *v;
    long i, j, ndim;
    long count;
    long nW, mW;
    double n, m, maze_size,core_iterations;

    double *store_y;
    long test_cell_num, test_source, test_target;
    double *test_dw;
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
	W  = (double* ) mxGetPr(prhs[0]);
    mW = mxGetM(prhs[0]);
    nW = mxGetN(prhs[0]);

	ww = (double* ) mxGetPr(prhs[1]);
	x  = (double* ) mxGetPr(prhs[2]);

    n = (long) mxGetScalar(prhs[3]);
    m = (long) mxGetScalar(prhs[4]);
    maze_size = (long) mxGetScalar(prhs[5]);
    core_iterations = (long) mxGetScalar(prhs[6]);

    /*additional agruments for testing derivatives, only change one cell's weight*/
    test_cell_num = (long) mxGetScalar(prhs[7]);
    test_source = (long) mxGetScalar(prhs[8]);
    test_target = (long) mxGetScalar(prhs[9]);
    test_dw = (double*) mxGetPr(prhs[10]);

    /*create outputs*/
	plhs[0] = mxCreateDoubleMatrix(1, (n+m)*maze_size*maze_size, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(n*maze_size*maze_size, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((n+m)*maze_size*maze_size, core_iterations,mxREAL);


	y = mxGetPr(plhs[0]);
    v = mxGetPr(plhs[1]);
    store_y = mxGetPr(plhs[2]);

    //mexPrintf("maze_size %f\n", maze_size);
    //mexPrintf("recurrent nodes n %f\n", n);
    //mexPrintf("input nodes m %f\n", m);
    //mexPrintf("core_iterations  %f\n", core_iterations);

    Maze_Forward_Calc(/*inputs*/W, mW, nW, ww, x,n, m, maze_size, core_iterations,
                       /*outputs*/ y, v, store_y, /*test args*/test_cell_num, test_source, test_target, test_dw);
}
