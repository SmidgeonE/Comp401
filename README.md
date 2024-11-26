# Boid Simulation

This is a C++ project using OpenMP and MPI to perform a multi-threaded N-Body simulation, with the choice of using it on a supercomputer of interconnected compute nodes.

# Requirements

You will need access to mpiicpx for compilation, as well as usual libs like OpenMP.

If some dependencies are missing, it should just automatically work if you do this on BC4 itself.

# Setup

First git clone it to your area. The actual simulation source and header files are in the programFiles folder. 

```
git clone https://github.com/SmidgeonE/Comp401.git

cd Comp401
```

I have provided bash scripts to make compilation easy. It is under the 'scripts' dir:

```
./scripts/compile.sh -h
```

I have used std::array among other things so details of the sim have to supplied at compile time with flags. Not all have to be supplied at compile time, but I didn't want to have to deal with compile time + run time arguments.
Here is an example of compiling an executable of 15 Boids with 4 OpenMP threads, and a simulation bounding box of size 30, that will save the simulation to a local .csv called BoidSimulation.csv:

```
./scripts/compile.sh -t 4 -n 15 -b 30 -o
```

This will produce a binary in the current directory with the name <numBoids>.exe . This can then be ran:

```
./15.exe
```

## Note

You can also add the separation, cohesion, and alignment constants at runtime as well to override the default values:

```
./15.exe <separationConstas> <cohesionConstant> <alignmentConstant>
```


Or, with MPI:

```
mpiexec -n <numMPIcores> ./15.exe
```

Alternatively, this can be ran as a job on BC4. There is an example BC4 script in ./bc4/ that you can use.

# Visualisation

You may need to install polars for the visualisation part: 

```
pip install polars
```


Once the BoidSimulation.csv file has been made in the base directory, you can create a gif that shows whats happening. From the base directory you can run:

```
python ./utils/visualise.py
```

# Examples:

Some premade gifs of the simulations are in ./visualisationExamples

It will take like 10-20 secs. Once its done, a simulation.gif file will appear in the base directory.


The default parameters will showcase a fairly normal behaviour: a sort of bifurcating flocking pattern.

![Simulation Example](./visualisationExamples/standard20.gif)

If you increase the separation parameter (by around 100x), you get this sort of spiralling flock behaviour.


![Simulation Example](./visualisationExamples/highSeparation.gif)



If you increase cohesion by about 10x, you get a very aggressive close flocking, which seems to in radius over time.

![Simulation Example](./visualisationExamples/highCohesion.gif)


High Alignment (3x default), will give a similar pattern, but it often decays to a straight line path for the flock. Low alignment will cause them to have a much more sporadic path.

![Simulation Example](./visualisationExamples/lowAlignment.gif)







