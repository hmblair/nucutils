## Overview

So far, there are four functions available:
   1. `getchar`, which gets the `n`-th character of the string, or all characters between `m` and `n` if two integers are passed;
   2. `getlength`, which gets the length of the string; and
   3. `embed`, which takes a FASTA file and converts it to an HDF5 file by   converting the bases into integers.
   4. `convert`, which takes a FASTA file and replaces all copies of the first specified base with the second specified base *in place*.

## Installation

To build, you will need installed copies of
   1. `samtools` and `HTSlib`.
   2. A copy of `HDF5` with parallel support enabled.
   3. `OpenMP` and `OpenMPI`.

These must be visible to `cmake`. If you are running on Sherlock, then 
```
ml load hdf5/1.14.4
ml load biology samtools/1.16.1
ml load cmake/3.24.2
```
is sufficient.

Installation is done by cloning the repository and running the installation script.
```
git clone https://github.com/hmblair/nucutils
cd nucutils
./configure
```
Don't forget to add the `./bin` directory to your path.


