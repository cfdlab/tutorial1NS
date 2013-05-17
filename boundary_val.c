#include "boundary_val.h"

/* ----------------------------------------------------------------------- */
/*                             Set Boundary Conditions                     */
/* ----------------------------------------------------------------------- */

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
                    int imax,
                    int jmax,
                    double **U,
                    double **V
                    ) {
    
    int i,j;
    /*Set boundary values along the columns*/
    for (j = 1; j <= jmax; j++){
        /*U velocities on right and left boundaries*/
        U[0][j] = 0;
        U[imax][j] = 0;
        /*V velocities on right and left boundaries (interpolated value)*/
        V[0][j]=-1*V[1][j];
        V[imax+1][j]=-1*V[imax][j];
    }
    
    /*Set boundary values along the rows*/
    for (i = 1; i <= imax; i++){
        /*V velocities on top and bottom boundaries*/
        V[i][0] = 0;
        V[i][jmax] = 0;
        /*U velocities on top and bottom boundaries (interpolated value)*/
        U[i][0]=-1*U[i][1];
        U[i][jmax+1]=2.0-1*U[i][jmax];
    }
    
}
