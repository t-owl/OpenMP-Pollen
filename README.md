# Overview

This Project requires knowledge of hared memory programming to parallelise and extend a given serial code You will recognise the basis of this project as being similar to that for [MPI-Pollen](https://github.com/t-owl/MPI-Pollen) but note there are some changes to initialisation and constants used in the simulation.

This code is a naïve model of pollen blowing in the wind and slowly falling to the ground. We initialise NUMPOLLEN pollen particles to a set of (x,y,z) and for each of MAXTIME (=3000) timesteps we mimic a wind function to give velocity in each axis which is used to update the (x,y,z) components of position. Finally, in that timestep we count the number of total pollen particles on the ground timestep. Updates on the amount of pollen on the ground is output every 50 timesteps.

The main aim of this project is to implement an efficient OpenMP code that models the falling pollen and performs the relevant best linear fit, while having a good performance.



## Parallel Code

To star with the parallelisation of the code the following steps were followed:

First it was a matter of identifying the parallel and serial regions of the code; the independent iterations, the independent variables, etc. Next was to set up tests, by splitting the program into different functions as well as reducing the problem size, this help to have a better understanding of the code and also to create a more efficient code. Finally once the programs worked correctly independently  it was a matter of putting everything together and scaling up the problem size. 

The parallelisation of the code is divided into two main parts, one part handles the simulation of pollen falling into the ground (Part 1), where as the other handles the computation of the best linear fit line (Part 2). 

**Part 1**. This part consists of 3 for loops, the first loop:

`for (timestep; timestep<MAXTIME; timestep++) {}`

It cannot be parallelised, because the results obtained for each time-step are dependent from the results obtain in the previous iteration. 

The second loop:

`for(i=0; i<LINEARARRAY; i++) {}`

It has many time consuming tasks such as `sqrt`, etc. In order to obtain a performance boost a `#pragma omp parallel for schedule(dynamic, 1024)` was added, this would schedule a dynamic assignment of iterations for each thread. This schedule is the preferred as some iterations would take more than than others. Specifically, the iterations where z[i] > 0.0 would take noticeably more time to compute than the iterations where this condition is not met (for this last case, no operations are performed). This kind of scheduling is appropriate when the iterations require different computational costs in order to balance the processing power between the set of threads. The fraction of iterations where no computation would be made (z[i] <= 0.0) would increase as time step increases, as more pollen would be sitting on the ground. A chunk size of 1024 is used, as this represents a good result obtained for timing measurements made.

Finally the third loop does not take as much time as the second one therefore the parallelisation of this loop is not necessary.

**Part 2.** On the other hand this part only consists of 1 loop:

`void linear_fit (double x_time[], double y_pollenGround[], double* m, double* c) {}`. 

A parallelisation scheme is taken on this function as well, taking care of the reduction clause which is used for accumulating variables. In this way, no racing condition occurs during the execution as OpenMP manages to handle these special operations for variables specified under the reduction parentheses.

It must be stated that this function is only executed once and its parallelisation has nearly no effect at the total processing time of this application. This has been modified for running in a parallel computing architecture for learning purposes.

## Set up & Run 

In this experiment two codes were produced, the serial solution and the parallel solution. 

### Serial solution

#### Set up

To start with the serial solution was a modified version of the original `pollen-v2.c` file, the modified version added one function two calculate the Y = mX+c formula. The Intel compiler module was needed as this was a simple c file. In order to compile and run the file a batch script was created to communicate with the `Slurm` system of `barkla`, this file contained the commands for **loading the Intel** module, and because the aim of this experiment is to explore and discuss parallelism all compilations were done without the compiler optimisation, therefore the **`-O0` flag** was used, finally to time the time taken by the program, the omp_get_wtime() function was used  inside the c program therefore **`-qopenmp flag`** was used along with **`OMP_NUM_THREADS=1`** this was implemented to set the number of threads used to 1.

```bash
## intel compiler
module load compilers/intel/2019u5 

# INTEL no-opt
echo INTEL no-opt
icc -qopenmp -O0 pollen-v2.c -o serial.exe
export OMP_NUM_THREADS=1; ./serial.exe
echo '-------'
```

#### Run

To run the bash file all it was needed in this case was the following command

```
sbatch run-serial.sh 
```

This would run the script described in the set up section, and create a `slurm` file as an output.

### Parallel solution

#### Set up

Similarly to the serial set up the **Intel c compiler** was needed, therefore it was loaded at the start of the created batch file, once it was loaded an if statement took care of checking if OpenMP was running in a single node, if this was the case the script would exit and terminate the program. After this checking the program was compiled using the **`O0` flag** as with the serial solution to further discuss parallelism, the **`-qopenmp`** flag was also declared to let the compiler that this code had OpenMP code.  Finally if the number of processes is not specified the program will run in 1 process, to specify the number of process the -c flag needs to be called along with the number of processes e.g. -c 4. This will be further explain in the run section of parallel solution.

```bash
# load modules
## intel compiler
module load compilers/intel/2019u5 


# check expected inputs (OpenMP is only supported on a single node)
if [ "$SLURM_NNODES" -gt "1" ]; then 
    echo more than 1 node not allowed
    exit
fi

# parallel using OpenMP
SRC=pollen-v2.c
EXE=${SRC%%.c}.exe
echo compiling $SRC to $EXE
icc -qopenmp -O0 $SRC -o $EXE && \
      (
# set number of threads
      export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} # if '-c' not used then default to 1
      echo using ${OMP_NUM_THREADS} OpenMP threads
      ./${EXE};echo
      ) \
      || echo $SRC did not built to $EXE
      


```

#### Run

In particular the OpenMP script needed the definition of the -c flag, this was done to declare the number of processes used, and this would determine on how the program was divided. to illustrate this two examples are presented. 

| Commands                   | Description                                                  |
| -------------------------- | ------------------------------------------------------------ |
| sbatch -c 2 run-openmp.sh  | to start a batch job while parsing "2" as the number of threads |
| sbatch -c 32 run-openmp.sh | to start a batch job while parsing "32" as the number of threads |

## Best Linear fit 

In order to find the best linear fit for the data produce by the simulation, the formula given in the *overview document* was used, where x is the time and y is the number of pollen that has fallen into the ground.

The line that determines the best linear fit is **Y = mX+c**. (Mathsisfun.com, 2017)

So to find m and c we apply the following formulas:

<img src="https://render.githubusercontent.com/render/math?math=m=\frac{N\sum_{xy}-\sum_{x}\sum_{y}}{N\sum_{x^{2}}-(\sum_{x})^{2}}#gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math=\color{White}m=\frac{N\sum_{xy}-\sum_{x}\sum_{y}}{N\sum_{x^{2}}-(\sum_{x})^{2}}#gh-dark-mode-only">

<img src="https://render.githubusercontent.com/render/math?math=c=\frac{\sum_{y}-m\sum_{x}}{N}#gh-light-mode-only"><img src="https://render.githubusercontent.com/render/math?math=\color{White}c=\frac{\sum_{y}-m\sum_{x}}{N}#gh-dark-mode-only">

In the implementation of the serial code one void function was created to calculate moth `m` and `c`, once the calculation was done m and c were returned to the main function to be printed into the terminal.

The code takes as input two arrays; one holds the total number of pollen that has fallen into the ground after each time-step as well, the other one holds the time-steps, both arrays take this data from the last 1000 steps. With this two inputs the functions are able to find m and c, and it is finally outputted in the console.

The resulting formula for this experiment was: 

y = 110.405819 x + 1995850.411865



To check the correctness of this formula a graph with the pollen data points was plotted against the formula in excel **Y = mX+c**.

![](https://i.imgur.com/m32VRZT.png)

![](https://i.imgur.com/6wJbjX9.png)

The first graph presents all the data taken from the experiment, the data shows an exponential growth of the number of pollen that has fallen into the ground, this growth happens from the time step 0 to the time step 1000 approx. After this a linear growth is present, more noticeable in the last 1000 steps

The second graph focuses on the last 1000 steps by using the data from the experiment a linear fit formula is estimated, which consequently as expected, the line of best fit goes along with the given data, and the estimated formula is close to the one calculated by the program, confirming that the derived formula is correct.

## Performance 

The performance of a parallel code is one of the priorities for this assignment, a good performance will show a good implementation of OpenMP functions.

The key aspects of performance are the speedup, efficiency and scalability.

Speed up is the ratio between the time taken in 1 core (or serial implementation) divided by the time taken in several cores (p)cores. i.e. the ideal speedup is going to be when running the program on the double number of cores, and the time takes half as long, so the speed up would be 2

Efficiency is the percentage of the speedup dived by the number of cores, so the ideal efficiency  would be 100%. Efficiency above 80% is considered good efficiency 

In this experiment the number of cores used to evaluate the performance were the following. Number of cores {1,2,3,4,6,8,10,13,16,19,32,40}.(Refer to appendix B for full performance data)

In the graph below we can see that the addition of threads does not affect the performance after the 15th thread

![](https://i.imgur.com/p0xnvpe.png)



In the graph below we can appreciate the speed up is not optimal, an orange line has been drawn to represent the ideal speed up as we can see it drops speedup drastically quite fast this might be due to not a great optimisation carried.

![](https://i.imgur.com/iMjpSIw.png)

Again in the graph below the efficiency graph goes as low as 12%

![](https://i.imgur.com/vux2v9w.png)

By looking at the graph we can conclude that the optimisation of the parallel program is not great, this could be directly linked to parts of the code not being parallel when they could, for instance one of this parts could be the first loop in (Part 1), In a previous assignment this part was parallel, due to a different requirement, in particular adding the best linear fit for the **last 1000 steps** and the way the code was structured, meant that parallelisation was not possible. Another reason could have been a poor implementation of the OpenMP code.

When looking at the data of the graphs we can see that the only good efficient (meaning efficiency above 80%) run was running with two threads. The quickest run is with 32 cores which took 37.663 seconds, but we are down to 16% of efficiency, which is not acceptable, therefore if we want to make a good use of our parallel program and reduce time we might use 2 cores instead which is the only run that has good efficiency. 

By looking at this data we can conclude the program does not follows a **strong scaling** approach as the number of cores increases the time to the solution does not decrease proportionally.

To get the time measured `omp_get_wtime()` was used through the code. variable called {t1,t2,t3,t4,t5} were used to time different parts of the code and by subtracting this values the running time of each section was found. 

**Accuracy**. Through all experiments the result for the **Y = mX+c** formula has been constant, `double` variables has been used across the code this means that is accurate to the 16 decimal place.

## MPI vs OpenMP

During the implementation of both programs, I have found OpenMP to be simpler and more readable almost as good as the serial implementation. Programming in OpenMP is relatively easy it just involves adding `pragma` directives, and the set number of threads the user wishes to run it on. 

MPI has overheads associated with the transferring of messages from one process to another, where as OpenMP has none of these overheads as threads can share variables. On the other hand OpenMP suffers from data racing (race condition), there are ways to mitigate this problem, but it can be troublesome.

MPI targets both distributed and shared memory systems, whereas OpenMP only targets shared memory systems

To run MPI programs a module made specifically for MPI has to be load, where as for OpenMP there is only need of the Intel/gcc compiler.

Finally the performance on MPI  has clearly been better through the assignments, this could be due to having spent more time working on MPI and therefore resulting in a better implementation, It is hard to judge completely as the approaches followed in each assignment  had small variations, still I would argue that MPI is more efficient in performance from the results obtained. 

In conclusion, in future implementations of parallel programming I would consider all of the points made, and for instance if an implementation involve distributed memory, then I would be forced to go with MPI, Otherwise I would have to consider all the points above.
