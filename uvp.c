#include "uvp.h"
#include <math.h>
#include <stdlib.h>

/* ----------------------------------------------------------------------- */
/*                             Function calculate_fg                       */
/* ----------------------------------------------------------------------- */

/**
 * Determines the value of F and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
 \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
 \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
 \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
 \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
 \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
 y^2} \right]_{i,j} \right) - \left[ \frac{\partial
 (uv)}{\partial x} \right]_{i,j} - \left[
 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
 \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */
void calculate_fg(
                  double Re,
                  double GX,
                  double GY,
                  double alpha,
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G
                  )

{
    
    int i ;
    int j ;
    double du2dx ;
    double duvdy ;
    double d2udx2 ;
    double d2udy2 ;
    
    double dv2dy ;
    double duvdx ;
    double d2vdx2 ;
    double d2vdy2 ;
        
    /*Determines the value of F according to the formula above with the help of temporary variables*/
    for ( i = 1 ; i < imax ; i++ )
    {
        for( j = 1 ; j <= jmax ; j ++ )
        {
            
            d2udx2 = ( U[i+1][j]  - 2*U[i][j] + U[i-1][j] ) / ( dx * dx) ;
            
            d2udy2 = ( U[i][j+1]  - 2*U[i][j] + U[i][j-1]) / (dy * dy )  ;
            
            du2dx = (1/dx) * ( ( (U[i][j] + U[i+1][j])/2 )*( (U[i][j] + U[i+1][j])/2 ) - ( (U[i-1][j] + U[i][j])/2 )*( (U[i-1][j] + U[i][j])/2 ) ) +
            alpha/dx * ( abs( U[i][j] + U[i+1][j] ) / 2  * ( U[i][j] - U[i+1][j] ) / 2 - abs( U[i-1][j] + U[i][j] ) / 2  * ( U[i-1][j] - U[i][j] ) / 2   ) ;
            
            duvdy = (1/dy) * ( ( V[i][j] + V[i+1][j] ) /2  *  ( U[i][j] + U[i][j+1] )/2 - (V[i][j-1] + V[i+1][j-1])/2 * (U[i][j-1] + U[i][j])/2  ) +
            alpha/dy * (abs( V[i][j] + V[i+1][j] ) /2  *  ( U[i][j] - U[i][j+1] )/2 - abs(V[i][j-1] + V[i+1][j-1])/2 * (U[i][j-1] - U[i][j])/2 ) ;
            
            
            
            F[i][j] = U[i][j]  + dt * ( 1/Re * ( (d2udx2 ) + (d2udy2) ) - (du2dx)  - duvdy + GX ) ;
            
            
        }
    }
    
    /*Determines the value of G according to the formula above with the help of temporary variables*/
    for ( i = 1 ; i <= imax ; i++ )
    {
        for( j = 1 ; j < jmax ; j ++ )
        {
            
            d2vdx2 = ( V[i+1][j]  - 2*V[i][j] + V[i-1][j] ) / ( dx * dx) ;
            
            d2vdy2 = ( V[i][j+1]  - 2*V[i][j] + V[i][j-1]) / (dy * dy )  ;
            
            
            duvdx = (1/dx) * ( ( V[i][j] + V[i+1][j] ) /2  *  ( U[i][j] + U[i][j+1] )/2 - (U[i-1][j] + U[i-1][j+1])/2 * (V[i-1][j] + V[i][j])/2  ) +
            alpha/dx * (( V[i][j] - V[i+1][j] ) /2  *  abs( U[i][j] + U[i][j+1] )/2 - abs(U[i-1][j] + U[i-1][j+1])/2 * (V[i-1][j] - V[i][j])/2 ) ;
            
            dv2dy = (1/dy) * ( ( (V[i][j] + V[i][j+1])/2 )*( (V[i][j] + V[i][j+1])/2 ) - ( (V[i][j-1] + V[i][j])/2 )* (V[i][j-1] + V[i][j])/2 )  +
            alpha/dy * ( abs( V[i][j] + V[i][j+1] ) / 2  * ( V[i][j] - V[i][j+1] ) / 2 - abs( V[i][j-1] + V[i][j] ) / 2  * (  V[i][j-1] - V[i][j]  ) / 2   ) ;
            
            
            
            G[i][j] = V[i][j]  + dt * ( 1/Re * ( (d2vdx2 ) + (d2vdy2) ) - (duvdx)  - dv2dy + GY ) ;
            
        }
    }
    
    /*Set boundary values along the columns*/
    for (j = 1; j <= jmax; j++){
        /*F values on right and left boundaries*/
        F[0][j] = U[0][j];
        F[imax][j] = U[imax][j];
    }
    
    /*Set boundary values along the rows*/
    for (i = 1; i <= imax; i++){
        /*G values on top and bottom boundaries*/
        G[i][0] = V[i][0];
        G[i][jmax] = V[i][jmax];
    }
    
}

/* ----------------------------------------------------------------------- */
/*                             Function calculate_rs                       */
/* ----------------------------------------------------------------------- */

/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **F,
                  double **G,
                  double **RS
                  ) {
	int i, j;
    for(i = 1; i <= imax; i++) {
        for(j = 1; j <= jmax; j++) {
            RS[i][j] = 1 / dt*( (F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy);
        }
    }
}

/* ----------------------------------------------------------------------- */
/*                             Function calculate_dt                       */
/* ----------------------------------------------------------------------- */

/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
                  double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V
                  ) {
    /*calculates maximum absolute velocities in x and y direction*/
    double umax=0, vmax=0;
    double a,b,c;
    int i, j;
    for(i = 1; i <= imax; i++) {
        for(j = 1; j<=jmax; j++) {
            if(abs(U[i][j])>umax)
                umax = abs(U[i][j]);
            
            if(abs(V[i][j])>vmax)
                vmax = abs(V[i][j]);
            
        }
    }
    
    /*Determines the minimum of dt according to stability criteria and multiply it by safety factor tau if tau is positive, otherwise uses the default value of dt*/
    if (tau>0){
        a = Re/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
        b = dx/umax;
        c = dy/vmax;
        if(a < b && a < c){
            *dt = tau * a;
        }
        else if(b < a && b < c){
            *dt = tau * b;
        }
        else{
            *dt = tau * c;
        }
    }
}

/* ----------------------------------------------------------------------- */
/*                             Function calculate_uv                       */
/* ----------------------------------------------------------------------- */

/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G,
                  double **P
                  ){
    int i;
    int j;
    /*Calculate the new velocity U according to the formula above*/
    for(i = 1; i <= (imax-1); i++){
        for(j = 1; j <= jmax; j++){
            U[i][j] = F[i][j]-(dt/dx)*(P[i+1][j]-P[i][j]);
        }
    }
    /*Calculate the new velocity V according to the formula above*/
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= (jmax-1); j++){
            V[i][j] = G[i][j]-(dt/dy)*(P[i][j+1]-P[i][j]);
        }
    }
}
