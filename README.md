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
