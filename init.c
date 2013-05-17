#include "helper.h"
#include "init.h"

/* ----------------------------------------------------------------------- */
/*                             Reading Parameters                          */
/* ----------------------------------------------------------------------- */


/**
 * This operation initializes all the local variables reading a configuration
 * file. For every variable a macro like READ_INT() is called passing it the
 * input filename and the variable to be written to. This macro calls
 * an operation read_int() augmenting the parameter set with the name of the
 * variable to be read. The read_int() operation parses the input file, extracts
 * the value of the variable, sets the variable and finally prints some debug
 * information. This is possible as the macro adds the name of the variable to
 * be set. All the helper operations can be found within helper.h and helper.c.
 *
 * @param szFileName char pointer to the filename
 * @param Re         Reynolds number
 * @param UI         initial velocity in  x-direction - used by init_uvp()
 * @param VI         initial velocity y-direction - used by init_upv()
 * @param PI         initial pressure - used by init_upv()
 * @param GX         gravitation x-direction
 * @param GY         gravitation y-direction
 * @param t_end      end time (not discrete in time steps)
 * @param xlength    domain length x-direction
 * @param ylength    domain lenght y-direction
 * @param dt         time step length: dividing t_end by dt gives the number of
 *                   time steps to perform. Actually dt is determined by a
 *                   function, so manipulating this value within the
 *                   configuration file should not affect the solution process
 *                   at all
 * @param dx         cell length x-direction
 * @param dy         cell length y-direction
 * @param imax       number of cells in x-direction
 * @param jmax       number of cells in Y-direction
 * @param alpha      uppwind-differencing-factor alpha
 * @param omg        relaxation factor omega
 * @param tau        safety parameter for time step calculation
 * @param itermax    max. number of pressure iterations
 * @param eps        tolerance limit for pressure calculation
 * @param dt_value   time steps for output (after how many time steps one should
 *                   write into the output file)
 */

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
                    /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
                    double *dt_value)           /* time for output */
{
    READ_DOUBLE( szFileName, *xlength );
    READ_DOUBLE( szFileName, *ylength );
    
    READ_DOUBLE( szFileName, *Re    );
    READ_DOUBLE( szFileName, *t_end );
    READ_DOUBLE( szFileName, *dt    );
    
    READ_INT   ( szFileName, *imax );
    READ_INT   ( szFileName, *jmax );
    
    READ_DOUBLE( szFileName, *omg   );
    READ_DOUBLE( szFileName, *eps   );
    READ_DOUBLE( szFileName, *tau   );
    READ_DOUBLE( szFileName, *alpha );
    
    READ_INT   ( szFileName, *itermax );
    READ_DOUBLE( szFileName, *dt_value );
    
    READ_DOUBLE( szFileName, *UI );
    READ_DOUBLE( szFileName, *VI );
    READ_DOUBLE( szFileName, *GX );
    READ_DOUBLE( szFileName, *GY );
    READ_DOUBLE( szFileName, *PI );
    
    *dx = *xlength / (double)(*imax);
    *dy = *ylength / (double)(*jmax);
    
    return 1;
}

/* ----------------------------------------------------------------------- */
/*                             Initializing U,V & P                        */
/* ----------------------------------------------------------------------- */

/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */

void init_uvp(
              double UI,
              double VI,
              double PI,
              int imax,
              int jmax,
              double **U,
              double **V,
              double **P
              )

{
    int i;
    int j;
    
    for ( i = 1 ; i<= imax;  i++ )
    {
        for (j = 1; j<= jmax; j++ )
        {
            
            U[i][j] = UI ;
            V[i][j] = VI ;
            P[i][j] = PI ;
            
        }
    }
    
}
