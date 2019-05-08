# cpd-particle-simulator
OpenMPI version of the code

## Copying file to the cluster (1st step)
    $ scp -R path/to/cpd-particle-simulator/mpi istTWXYZ@cluster.rnl.tecnico.ulisboa.pt:/mnt/cirrus/users/Y/Z/istTWXYZ

## Run
    $ ssh istTWXYZ@cluster.rnl.tecnico.ulisboa.pt
    $ Password:
    $ istTWXYZ@borg:~$
    $ make
    $ srun -n <number of nodes> simpar-mpi <args>
