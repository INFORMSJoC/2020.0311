# 2020.0311

[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# SPGP

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[A Simulation Optimization Approach for Appointment Scheduling Problem with Decision-Dependent Uncertainties](https://doi.org/) 
by T. Homem-de-Mello, Q. Kong and R. Godoy-Barba. 


## Cite

To cite this software, please cite the [paper](https://doi.org/) using its DOI.

Below is the BibTex for citing this version of the code.

```
@article{CacheTest,
  author =        {T. {Homem-de-Mello}, Q. Kong and R. {Godoy-Barba}},
  publisher =     {INFORMS Journal on Computing},
  title =         {A Simulation Optimization Approach for Appointment Scheduling Problem with Decision-Dependent Uncertainties},
  year =          {2022},
  doi =           {https://doi.org/},
  url =           {https://github.com/INFORMSJoC/2020.0311},
}  
```

## Description

The goal of this software is to solve the problem of appointment scheduling with time-dependent no-shows. In such problem the goal is to determine the 
appointment time of a list of patients, taking into account that each patient may not show up for their appointment. Moreover, the probability of no-show
depends on the time scheduled for the appointment. This is therefore an optimization problem with decision-dependent uncertainty.

## Running the code

The code runs in Matlab. The file `SPGP.m` contains a script with all the parameters of the problem, and then a call to a routine that implements the
optimization algorithm. The default values of the parameters can be changed, following the instructions on the script file.

There are also two parameters for the algorithm itself. One is the initial sample size used by the simulation-based method, and the other (`Trials`) is the desired number of runs. As described in the paper, a multi-start approach is implemented to deal with the non-convexity of the problem, so each run randomly selects a feasible starting point for the algorithm.

Once the code is run, the results (optimal schedule, and a 95% confidence interval for the optimal cost) are saved on files called `schedule_<type>_<run>.mat`, where `<type>` indicates the show-up probability function class (1, 2 or 3) and `<run>`is the run number.

## Results

Figure 3 in the paper, reproduced below, shows the optimal schedules for 6 different show-up probability functions, under exponentially distributed service times. In the script, one selects one of three classes of  show-up probability functions: linear, quadratic, or cosine. The parameter `P1` corresponds to the show-up probability at time 0, whereas the meaning of `P2` depends on which of the three classes of show-up probability functions is being used. 


![Figure 1](https://github.com/thmello1/2020.0311/files/8496907/schedule_p_esp.pdf)


## Replicating

To replicate approximately the results in the figure, run the `SPGP.m` script changing the value of `ShapeShowUpProb` to the corresponding class, and the values of `P1` and `P2` accordingly. The plots can be obtained using the script `plot_schedule.m`. Note that it is not possible to replicate exactly the results in the figure since the algorithm uses random samples and random starting points. Still, the results should looks similar to those in the figure.

