#include <stdio.h>
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
double xDelta, yDelta, r;                           // wind changes
double eps=1.0E-06;
int numGround;

int initialise();

void linear_fit (double x_time[], double y_pollenGround[], double* m, double* c){
    
    int i;
    double sumXY=0;
    double sumX=0;
    double sumX2=0;
    double sumY=0;
    for(i=0;i<LINEARARRAY;i++){
        sumXY=sumXY+x_time[i]*y_pollenGround[i];
        sumX=sumX+x_time[i];
        sumY=sumY+y_pollenGround[i];
        sumX2=sumX2+x_time[i]*x_time[i];
    }
    sumXY=sumXY/LINEARARRAY;
    sumX=sumX/LINEARARRAY;
    sumY=sumY/LINEARARRAY;
    sumX2=sumX2/LINEARARRAY;

    *m=(sumXY-sumX*sumY)/(sumX2-sumX*sumX);
    *c=(sumX2*sumY-sumXY*sumX)/(sumX2-sumX*sumX);
}


int main(void) {
  int rootRank=0, timestep=0;
  int i;
  int numThreads=-1;
  int num=NUMPOLLEN;
  
  int counter=0;
  double x_time[LINEARARRAY];
  double y_pollenGround[LINEARARRAY];

  printf("Initialising...");
  int rc = initialise();  //rc=0 good, rc<0 bad
  if (rc<0) {
    exit(rc);
  }
  printf("   DONE\n");

  double t0 = omp_get_wtime();

  for (timestep; timestep<MAXTIME; timestep++) {
      for (i=0; i<num; i++) {
	if(z[i] > 0.0) {
	  /*
	    change in component of velocity depends on position from eye of storm
	  */
	  xDelta = abs(eye_x - x[i]);
	  yDelta = abs(eye_y - y[i]);
	  r = sqrt(xDelta*xDelta + yDelta*yDelta)/200.0;  // a form of normalisation
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
      for (i=0; i<num; i++) {
	if(z[i] <= 0.0) numGround++;
      }

      if (timestep%50==0 || timestep==MAXTIME-1) printf("Timestep %d: %d particles on ground\n", timestep, numGround);

    if (timestep >=2000 && timestep <=2999 ){
        x_time[counter]= timestep;
        y_pollenGround[counter] = numGround;

        counter += 1;
    }

  } // end time

  /* 
     Note that you should include within the total time, 
     that taken to find best fit for last 1000 points 
  */


  double m;
  double c;

  double secondsTaken = omp_get_wtime() - t0;
  printf("%d pollen for %d timesteps on %d threads takes %f seconds\n", NUMPOLLEN, MAXTIME, numThreads, secondsTaken);
  
  linear_fit(x_time,y_pollenGround, &m, &c);
  printf("y = %lf x + %lf\n\n",m,c);

  } // main

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
      eye_x = 0.0; eye_y=0.0;


      return 0; // good
    } 
}
