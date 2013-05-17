#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>

/* CFD Lab - Worksheet 1 - Group 3
 * Camacho Barranco, Roberto
 * Gavranovic, Stefan
 * Valizadeh, Mahyar
 */

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){
    /*Geometry data*/
    double xlength;		 /*domain size in x-direction*/
    double ylength;		/*domain size in y-direction*/
    /* The computational domain is thus = [0; xlaenge] X [0; ylaenge].*/
    int imax;		/* number of interior cells in x-direction*/
    int jmax;		/* number of interior cells in y-direction*/
    double dx;		/* length dx of one cell in x-direction*/
    double dy;		/* length dy of one cell in y-direction*/
    /* Time-stepping data:*/
    double t;		/* current time value*/
    double t_end;		/* final time tend*/
    double dt;		/* time step size dt*/
    double tau;		/* safety factor for time step size control T*/
    double dt_value;	/* time interval for writing visualization data in a file*/
    int n;			/* current time iteration step*/
    /* Pressure iteration data:*/
    int itermax;		/* maximum number of pressure iterations in one time step*/
    int it;			/* SOR iteration counter*/
    double res;		/* residual norm of the pressure equation*/
    double eps;		/* accuracy criterion epsilon (tolerance) for pressure iteration (res < eps)*/
    double omg;		/* relaxation factor omega for SOR iteration*/
    double alpha;		/* upwind differencing factor alpha (see equation (4))*/
    /* Problem-dependent quantities:*/
    double Re;		/* Reynolds number Re*/
    double GX,GY;		/* external forces gx; gy, e.g. gravity*/
    double UI,VI,PI;	/* initial data for velocities and pressure*/
    /* Arrays*/
    double **U;		/* velocity in x-direction*/
    double **V;		/* velocity in y-direction*/
    double **P;		/* pressure*/
    double **RS;		/* right-hand side for pressure iteration*/
    double **F,**G;		/* F;G*/
    int n_div;
    
    /* read the program configuration file using read_parameters()*/
    read_parameters("cavity100.dat", &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value);
    
    /* set up the matrices (arrays) needed using the matrix() command*/
    U = matrix(0, imax+1, 0, jmax+1);
    V = matrix(0, imax+1, 0, jmax+1);
    P = matrix(0, imax+1, 0, jmax+1);
    RS = matrix(0, imax+1, 0, jmax+1);
    F = matrix(0, imax+1, 0, jmax+1);
    G = matrix(0, imax+1, 0, jmax+1);
    
    /* initialize current time and time step*/
    t = 0;
    n = 0;
    
    /* create the initial setup init_uvp()*/
    init_uvp(UI, VI, PI, imax, jmax, U, V, P);
    
    
    /* ----------------------------------------------------------------------- */
    /*                             Performing the main loop                    */
    /* ----------------------------------------------------------------------- */
    
    /* Upperbound t_end+dt/10 to be sure that it runs for t=t_end */
    while (t<t_end+dt/10){
        /*	Select dt*/
        calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
        /*	Set boundary values for u and v according to (14),(15)*/
        boundaryvalues(imax, jmax, U, V);
        /*	Compute F(n) and G(n) according to (9),(10),(17)*/
        calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V ,F , G);
        /*	Compute the right-hand side rs of the pressure equation (11)*/
        calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);
        /*	Set it := 0*/
        res = 1.0;
        it = 0;
        /*	While it < itmax and res > eps*/
        while(it < itermax && res > eps){
            /*	Perform a SOR iteration according to (18) using the*/
            /*	provided function and retrieve the residual res*/
            sor(omg, dx, dy, imax, jmax, P, RS, &res);
            /*	it := it + 1*/
            it++;
        }
        /*	Compute u(n+1) and v(n+1) according to (7),(8)*/
        calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);
        /*	Output of u; v; p values for visualization, if necessary*/
        
        n_div=(dt_value/dt);
        if(n % n_div == 0){
            write_vtkFile("cavity100", n , xlength, ylength, imax, jmax, dx, dy, U, V, P);
        }
        /*	t := t + dt*/
        t = t + dt;
        /*	n := n + 1*/
        n++;
        
    }
    
    
    /* Destroy memory allocated*/
    free_matrix(U, 0, imax+1, 0, jmax+1);
    free_matrix(V, 0, imax+1, 0, jmax+1);
    free_matrix(P, 0, imax+1, 0, jmax+1);
    free_matrix(RS, 0, imax+1, 0, jmax+1);
    free_matrix(F, 0, imax+1, 0, jmax+1);
    free_matrix(G, 0, imax+1, 0, jmax+1);
    return -1;
}
