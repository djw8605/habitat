Habitat Simulation
==================


Habitat is a parallel multi-agent simulator.


Building
--------

Building is done with the makefile.  You can start the build with the command:

    $ make all

By default, it will make the executable habitat.exe.


Running
-------

Two configuraiton files control how habitat is executed, habitat.in and 
makeland.in.  Both files must be available in the current directory when
starting habitat.exe.

Parallelization of habitat.exe is done with OpenMP.  By default, habitat.exe
will run as many threads as there are cores on the executing machine.  The
number of OpenMP threads can be controlled through the environment variable
`OMP_NUM_THREADS`.  Set it in the environment with:

    $ export OMP_NUM_THREADS=8

The above command will allow OpenMP to start 8 threads on the current machine.

Running of habitat.exe is by simplying starting it:

    $ ./habitat.exe



