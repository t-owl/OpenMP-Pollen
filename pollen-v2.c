#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdlib.h>  /* for rand() */
#include <math.h>

// NUMPOLLEN needs to be a square number (see initialisation)
#define NUMPOLLEN 4000000
#define FALLRATE 0.0005
#define MAXTIME 3000
#define STORMFORCE 0.0008
#define LINEARARRAY 1000

// global vars
double vx[NUMPOLLEN], vy[NUMPOLLEN], vz[NUMPOLLEN];
double x[NUMPOLLEN], y[NUMPOLLEN], z[NUMPOLLEN];
double eye_x, eye_y;                                // position of centre of storm
double eps=1.0E-06;
int numGround;

int initialise();

void linear_fit (double x_time[], double y_pollenGround[], double* m, double* c) {

    int i;
    double sumXY=0;
    double sumX=0;
    double sumX2=0;
    double sumY=0;

    // parallelise for loop and take care of the reduction which is used for accumulating variables.
    #pragma omp parallel for reduction(+:sumXY,sumX,sumY,sumX2)
    for(i=0; i<LINEARARRAY; i++) {
        sumXY+=x_time[i]*y_pollenGround[i];
        sumX+=x_time[i];
        sumY+=y_pollenGround[i];
        sumX2+=x_time[i]*x_time[i];
    }
    sumXY=sumXY/LINEARARRAY;
    sumX=sumX/LINEARARRAY;
    sumY=sumY/LINEARARRAY;
    sumX2=sumX2/LINEARARRAY;

    *m=(sumXY-sumX*sumY)/(sumX2-sumX*sumX);
    *c=(sumX2*sumY-sumXY*sumX)/(sumX2-sumX*sumX);
}


int main(void) {
    int i;
    int rootRank=0, timestep=0;
    int counter=0;
    double x_time[LINEARARRAY];
    double y_pollenGround[LINEARARRAY];

    int numThreads=atoi(getenv("OMP_NUM_THREADS"));
    /*
      The number of running threads for OpenMP sections of the
      code could be overwritten when using omp_set_num_threads()
      calls as well.
      By default, OpenMP uses the number of threads specified by the
      OMP_NUM_THREADS environment variable if this is set.
     */

    double t4 = omp_get_wtime();
    printf("Initialising...");
    int rc = initialise();  //rc=0 good, rc<0 bad
    if (rc<0) {
        exit(rc);
    }
    printf("   DONE\n");
    double t5 = omp_get_wtime();

    // Start measuring processing time
    double t1 = omp_get_wtime();
    //(First loop)
    for (timestep; timestep<MAXTIME; timestep++) {
     
        //(Second loop)Dynamic pragma for
        #pragma omp parallel for schedule(dynamic, 1024)
        for (i=0; i<NUMPOLLEN; i++) {
            if(z[i] > 0.0) {
                /*
                  change in component of velocity depends on position from eye of storm
                */
                double xDelta = abs(eye_x - x[i]);
                double yDelta = abs(eye_y - y[i]);
                double r = sqrt(xDelta*xDelta + yDelta*yDelta)/200.0;  // a form of normalisation
                vx[i] += STORMFORCE * r;
                vy[i] += STORMFORCE * r;
                vz[i] += -FALLRATE/(r*r+eps); // prevent divide-by-zero

                /* useful debugging output  */
                //printf("step %d particle %d at (%g,%g,%g), r=%g, now has velocity (%g,%g,%g)\n", \
                //timestep, i, x[i], y[i], z[i], r, vx[i], vy[i], vz[i]);

                x[i] += vx[i];
                y[i] += vy[i];
                z[i] += vz[i];
            }
        }

        // determine # on ground per process
        numGround=0;
        //(Third loop)
        for (i=0; i<NUMPOLLEN; i++) {
            if(z[i] <= 0.0) numGround++;
        }

        if (timestep%50==0 || timestep==MAXTIME-1) printf("Timestep %d: %d particles on ground\n", timestep, numGround);

        if (timestep >=2000 && timestep <=2999 ) {
            x_time[counter]= timestep;
            y_pollenGround[counter] = numGround;

            counter += 1;
        }

    }

    /*
       Note that you should include within the total time,
       that taken to find best fit for last 1000 points
    */
    

    double m;
    double c;

    double t2 = omp_get_wtime();
    linear_fit(x_time,y_pollenGround, &m, &c);
    double t3 = omp_get_wtime();
    
    printf("\n\ny = %lf x + %lf\n\n",m,c);

    /*
     * End measuring time, fitting time has been included as well
     * in order to comply with the comments below.
     */
    

    double secondsTaken = t3 - t1;
    printf("%d pollen for %d timesteps on %d threads takes %f seconds\n", NUMPOLLEN, MAXTIME, numThreads, secondsTaken);
    printf("Time taken to compute slop and intercept: %g secs\n", t3 - t2 );
    printf("Time taken to Initialse variables: %g secs\n", t5 - t4 );

}


int initialise() {
    /*
       domain is [-100,100] in x,y and [0,50] for z
       uniform distribution throughout (x,y,z) positions
       all initial velocities are zero
    */
    int i, numRows, numCols;
    numRows=sqrt(NUMPOLLEN);
    numCols = numRows;
    if (numRows*numCols != NUMPOLLEN) {
        printf("Error num needs be square - abort\n");
        return -1;
    }
    else {
        int row, col;
        double boxwidth = 200.0;
        double sep = boxwidth/(double)numRows;
        for (i=0; i<NUMPOLLEN; i++) {
            row = i/numRows;
            col = i - row*numRows;
            // printf("NUM=%d, i=%d, row,col=%d,%d\n", NUMPOLLEN, i, row, col);
            x[i] = -0.5*boxwidth + sep*(double)col;
            y[i] = -0.5*boxwidth + sep*(double)row;
            z[i] =  240.0;
            vx[i] = 0.001;
            vy[i] = (i%2==0) ? 0.001 : -0.001;
            vz[i] = (i%10==0) ? 1.0  : 2.0;
        }

        // set centre of storm
        eye_x = 0.0;
        eye_y=0.0;


        return 0; // good
    }
}