# ParallelUncertaintyHyperthermia
## Evaluates the influence of uncertainties in the temperature distribution 
<hr/>
<p>This is a parallelized version of the <a href='https://github.com/antonioMarchese/SerialUncertaintyHipertermia'>SerialUncertaintyHyperthermia</a> repository code, in which it is possible to analyze the influence of uncertainties on experimental parameters in hyperthermia treatment. The parallelization strategy allows many experiments to be performed while keeping a low execution time, thus obtaining better results.</p>
<hr/>

## Parallelization
<p>The parallelization strategy consists of two steps:</p>
<ol>
  <li>Parallelize the execution of the experiments</li>
  <li>Parallelize the resolution of the differential equation</li>
</ol>

### First step
<p>In the first stage, MPI is used to parallelize the execution of the 'experiment loops'. MPI (Message Passing Interface) manages a parallel computation on a distributed memory system. Thus, several processes are used, each one of them responsible for a specific number of experiments. This allows the code to be 'dismembered' into several pieces, with tasks being performed concurrently. However, it is necessary to carry out communication between the different processes, which often leads to a decrease in code efficiency.</p>

### Second step
<p>OpenMP is used to parallelize the differential equation solving. OpenMP API provides a relaxed-consistency, shared-memory model. In this step, each processor divides the task of executing the 'space loops' into its threads. In this step, each processor divides the task of executing the 'space loops' into its threads. This requires extra attention because the tasks to be parallelized cannot depend on each other. </p>
<hr/>

## To do

<p>First of all, you need must create a 'results' folder in your own repository. Otherwise, it is possible that you will get a segmentation fault, as the program's results are written to this folder in order to maintain a certain organization.</p>
<p>Then, all you got to do is (make sure you have both MPI and OpenMP installed in your machine): </p>
<h3>Compile: </ h3>
<h4>mpiCC -lm -fopenmp -O3 main.cpp -o 'file_name' </h4>
<h3>Run: </ h3>
<h4>mpirun -np 'number_of_processes' ./'file_name' </h4>
