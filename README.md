## Atrophy_FA_integration ##

This script expands the atrophy code by using voxel-wise shear modulus values instead of region-wise shear moduli.
The Mu values are currently calculated with 
  mu_element = (-(fa_value / 0.0037) + 182.4)*1e-6; 

Input files:
- mesh .inp file
- FA .inp file
- .prm file

## Installation

The efi library requires deal.II 9.1.1 to be installed with trilinos enabled. It is recommended to install dealii via spack. Further, the boost libraries *filesystem*, *tti*, and *any* are required. 

1. Make project directory, e.g. *myproject*.
2. Enter the directory.
3. Copy the src folder to the project directory.
4. Create a build directory.
5. Enter the build directory.
6. Run cmake (out-of-source-build).

To do so, run the following comands in a terminal:
<pre><code>
spack load dealii@...
spack load boost@...
spack load openmpi@...
mkdir myproject
cd myproject
cp -a /path/to/src .
mkdir build
cd /build 
cmake ../src
make debug
</code></pre>
or alternatively
<pre><code>
make release
</code></pre>

## Runnig the executable
The executable can be run form the terminal via
<pre><code>
mpirun -np 1 --bind-to socket efi_vlab 3 myparams.prm
</code></pre>
