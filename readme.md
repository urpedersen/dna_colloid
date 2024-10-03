# Model of DNA functionalized nano-particles

This program simulates a model of DNA functionalized nano-particles.
The program is written in C++,  and developed on a Linux system. 

The model is a coarse-grained representation,
where a colloidal particles are represented as soft repulsive-spheres, and 
DNA strands are represented as "arms" that "hold hands" with one other
if an angle and distance criteria is meet. The DNA arms energies 
are modeled as truncated harmonic potentials.

Trajectories are generated using Monte Carlo simulations, where the
orientation of the particles represented using quaternions. Simulations can 
be done in the canonical ensemble (NVT), grand canonical ensemble (µVT), or in the
isothermal-isobaric ensemble (NpT). The latter can be done using different
types of volume changes such as isotropic, anisotropic, or only in z-direction.
Finally, simulations can be done with a harmonic potential bias potential on the
number of particles, N, kappa(N-N_0)^2.

Thermodynamic information is written to a file. 
The particle trajectories are saved in the `xyz` format and as VMD files.
A range of orientation order-parameters can be calculated on the fly.
These include the Legendre polynomials, qubatic, and Hexgonal-diamond.

## Model details

Let $N$ be the number of particles, $Z$ the number of strands per particle, 
$r_{ij}$ the vector pointing from particle $i$ to $j$. 
Let $\beta^{-1} = k_B T$. The Hamiltonian is expressed as:

$$
\beta U = \frac{1}{2} \sum_{i=1}^N \left( \sum_{j=1}^N u_{ij}^{\text{r}} + \sum_{a=1}^Z u_{ia}^{\text{b}} \right)
$$

Here $u_{ij}^{\text{r}}$ is repulsive steric interactions between colloids and strands, and
$u_{ia}^{\text{b}}$ is attractive bonding between strands.

### Steric repulsions

Steric interactions are modeled as repulsive inverse power laws:

$$
u_{ij}^{\text{r}} = \left( \frac{\sigma_\text{r}}{r_{ij}} \right)^n + S(r_{ij})
$$

where $n$ controls the hardness and $\sigma_\text{r}$ defines the range. 
The shift function $S(r_{ij})$ truncates the pair repulsion at $r_\text{c}$ as follows:

$$
S(r_{ij}) =
\begin{cases}
-\left( \frac{\sigma_\text{r}}{r_\text{c}} \right)^n & \text{if } r_{ij} < r_\text{c} \\
-\left( \frac{\sigma_\text{r}}{r_{ij}} \right)^n & \text{otherwise}
\end{cases}
$$

### Bonding Strands

Bonds between strands are modeled as truncated springs. 
Let $v_{ia}$ be a normal vector from the center of particle $i$ in the direction of the $a$-th strand, 
and $\sin^2(\theta_{ija}) = 1 - (r_{ij} \cdot v_{ia})^2$. The bond distance $d_{iajb}^2$ is given by:

$$
d_{iajb}^2 = \frac{(r_{ij} - \sigma_\text{b})^2}{r_\text{max}^2} + \frac{\sin^2(\theta_{ija}) + \sin^2(\theta_{jib})}{\sin^2(\theta_\text{max})}
$$

The strand bond energy is a truncated harmonic potential (truncated spring):

$$
u^\text{b}_{ia} = \varepsilon \sum_{j=1}^N \sum_{b=1}^Z (d_{iajb}^2 - 1) H(1 - d_{iajb}^2) h_{iajb}
$$

where $H(x)$ is the Heaviside function, and $h_{iajb} = 1$ if strand $ia$ bonds with strand $jb$, and $h_{iajb} = 0$ otherwise.

Strands may only bind once, $\sum_{jb} h_{iajb} = [0, 1]$, and bonding is consistent, $h_{iajb} = h_{jbia}$. The bonding matrix $h_{iajb}$ is chosen to minimize the total bonding energy.

## Build from source

### Get the source code

Clone the repository hosted at GitHub with

```sh
git clone https://github.com/urpedersen/dna_colloid.git
```

The source code files are located in the `dna_colloid` directory.

```sh
cd dna_colloid
```

### Build using make (linux)
Build with make (in directory `dna_colloid` with the `makefile`)

```sh
make
```

The executable will be in the `bin` directory.

Later, you can clean directories with `make clean`.

## Execute the program

### My first run
After the program is built, you should be able to run the examples
located in the `./examples` folder. 

```sh
cd examples/default
../../bin/dna_colloid
```

You should see the program output a lot of information about the simulation.

### Input files

#### Variables
The program reads input from a file or the command line 
(with priority to the input given in the command line).
The input file is a text file with a list of parameters and values,
see the `examples/default/variables.txt` file for an example.
Here is how the first lines could look:

```
--verbose=9
--seed=1
--kT=1.0
...
```

#### Arms
The program also reads the DNA arms from a file named `arms.xyz`,
in the `xyz` format.
As an example, a tetrahedral DNA nano-particle with 4 arms could look like this:

```
4
Tetrahedron
0   1.0  1.0  1.0
0  -1.0 -1.0  1.0
0  -1.0  1.0 -1.0
0   1.0 -1.0 -1.0
```

More examples can be found in the `arms` directory.

#### Input configuration

An input configuration of particles position and orientation 
can be given in the `ini.xyz` file in the `xyz` format.
The first lines could look like this:

```
6
MC simulation: time=0 bboxX=5 bboxY=5 bboxZ=5 beta=1e+08 kT=1e-08 V=125 columns={type,x,y,z,ix,iy,iz,q0,q1,q2,q3,u}
0 2.22183 4.43218 4.55494 0 0 0 0.485023 0.503605 0.495247 0.515621 -5.85583
0 4.12581 1.03037 3.92471 0 0 0 0.491939 0.508109 0.494759 0.50501 -6.79104
0 4.14594 2.12176 0.175377 0 0 1 0.491279 0.513913 0.500302 0.494202 -5.70302
0 1.61161 4.98029 0.168304 0 0 0 0.520758 0.441036 0.521532 0.512155 -4.80811
0 2.25278 4.43165 2.05167 0 0 0 0.514068 0.487572 0.50493 0.493005 -5.86147
0 3.48872 4.43875 3.30815 0 0 0 0.515072 0.477444 0.494207 0.512355 -6.17234
```

The first lines contain the number of particles.
The second line contains a comment. Here,
`time` is the simulation time, `bboxX`, `bboxY`, `bboxZ` is the box size,
`beta` is the inverse temperature, `kT` is the thermal energy,
`V` is the volume, and `columns` is the column names described below.

Each following line contain one particle with the columns
`type,x,y,z,ix,iy,iz,q0,q1,q2,q3,u`.
Here, where `x`, `y`, `z` is the positions,
`ix`, `iy`, `iz` is the image coordinates,
`q0`, `q1`, `q2`, `q3` is the orientation quaternion,
and `u` is the single-particle potential-energy.

### Output

#### The pair potential
The pair potential is written to the files `potential_pair.dat`,
`potential_arm.dat`, `potential_arm_angle.dat` and 
`potential_arm_distance_angle.dat`. 
Output can be tweaked with the `print_pair_num_points` and
`print_angle_num_points` input variables.

#### Thermodynamic output

The program writes the thermodynamic output to the `ener.dat` file.
Use the `print_ener_freq` input variable to set the frequency of the output.
The columns in the ascii file are

1. Time (steps)
2. Potential (internal) energy
3. Number of particles
4. Volume
5. Box Length in x-direction
6. Box Length in y-direction
7. Box Length in z-direction
8. Orientational-order parameter (Q)
9. Orientational-order bias energy
10. Particle bias energy 
11. Total energy

#### Configuration output

The program writes the final configuration to the `end.xyz`.
Use the `print_conf_freq` input variable to set the frequency of the output.
This can be used to restart the simulation, by moving the `end.xyz` to
and `ini.xyz` in a different directory. The `conf.xyz` is a trajectory of positions and orientations.

##### VMD output

The program writes the final configuration to the `end.vmd` file, and
trajectories as `conf0.vmd, conf1.vmd, ...` files. These can be visualized with VMD.
Use the `print_vmd_freq` input variable to set the frequency of the output.

##### Network

The output file `network.dat` contains the network of connected particles.
Below is an example of the file format:

```
54 160 100
# Above: [ number of connections ; number of nodes ; time ]. Below: [ particle ; arm ; energy ; connecting particle ; connecting arm ; is connected ]
3 0 -0.364344 72 1 1
5 0 -0.0862437 17 1 1
6 0 -0.799204 80 1 1
7 1 -0.110971 8 0 1
8 0 -0.110971 7 1 1
9 1 -0.854967 60 0 1
10 0 -0.267134 129 1 1
13 1 -0.098108 151 2 1
...
```

#### Standard output

The standard output of the default run is a lot of information about the simulation:

```
Running ../../bin/dna_colloid
Read variable:  --verbose=9
Read variable:  --seed=1
Read variable:  --kT=1 
Read variable:  --chem_pot=0 
Read variable:  --time=0 
Read variable:  --bboxX=6 
Read variable:  --bboxY=6 
Read variable:  --bboxZ=6 
Note: For Lennard-Jones potential, put [ pair_eps=0 ; arm_eps=0 ; pair_12_eps=4 ; pair_6_eps=-4 ] or [ pair_eps=0 ; arm_eps=0 ; pair_12_eps=1 ; pair_6_eps=-2 ].
Read variable:  --pair_eps=1 
Read variable:  --pair_sigma=1 
Read variable:  --arm_eps=1 
Read variable:  --arm_r_center=1 
Read variable:  --arm_r_width=0.2 
Read variable:  --max_angle=40 
Max angle: theta_max = 40 degree = 0.698132 radians.
Read variable:  --pair_12_eps=0 
Read variable:  --pair_6_eps=0 
Read variable:  --pair_4_eps=0 
Note: Pair potential: pair_eps*(pair_sigma/r)^18
potential.set_parameters: pair_eps = 1
potential.set_parameters: pair_sigma = 1
potential.set_parameters: arm_eps = 1
potential.set_parameters: arm_r_center = 1
potential.set_parameters: arm_r_width = 0.2
potential.set_parameters: arm_sin_width = 0.642788
potential.set_parameters: pair_12_eps = 0
potential.set_parameters: pair_6_eps = 0
potential.set_parameters: pair_4_eps = 0
Read variable:  --cut_distance=1.2 
Note: only_bind_to_lowest_arm_energy = {0,1,2} 
Read variable:  --only_bind_to_lowest_arm_energy=0
Note: Bind to all arms. 
Read variable:  --num_extra_sticky_spots=0
Do not add extra sticky spots.
Note:   Selections of order parameters:                                                       
Note:   type_of_orientational_order=0    ->    Do not calculate order parameter               
Note:   type_of_orientational_order=1    ->    Use legendre4                                  
Note:   type_of_orientational_order=2    ->    Use legendre4_tagged      orientation relative to particle 0 
Note:   type_of_orientational_order=3    ->    Use qubatic               as in [JPC-B (2011) 115,14205]       
Note:   type_of_orientational_order=4    ->    Use legendre4_max_bonds   with Nb=4           
Note:   type_of_orientational_order=6    ->    Use ...                   with Nb=6           
Note:   type_of_orientational_order=8    ->    Use ...                   with Nb=8           
Note:   type_of_orientational_order=14   ->    Use legendre4_fix_bonds   with Nb=4           
Note:   type_of_orientational_order=16   ->    Use ...                   with Nb=6           
Note:   type_of_orientational_order=18   ->    Use ...                   with Nb=8           
Note:   type_of_orientational_order=24   ->    Use legendre4_neighbours  with Nb=4           
Note:   type_of_orientational_order=26   ->    Use ...                   with Nb=6           
Note:   type_of_orientational_order=28   ->    Use ...                   with Nb=8           
Note:   type_of_orientational_order=106  ->    Use hexdia (optimized for hexagonal diamond)  
Read variable:  --type_of_orientational_order=1
Note: If orientational_order_threshold is enables (1) then the order is unity if it is within the threshold value and zero otherwise.
Read variable:  --enable_orientational_order_threshold=0
Read variable:  --orientational_order_threshold_max=1.25 
Read variable:  --orientational_order_threshold_min=0.75 
Read variable:  --umbrella_orientational_order_kappa=0 
Read variable:  --umbrella_orientational_order_center=0 
Note: Orientational order umbrella is disabled (kappa<=0). 
Read variable:  --skin_distance=0.3 
Read variable:  --num_added_particles=0
Read variable:  --print_random_numbers=0
Read variable:  --print_pair_num_points=250
Read variable:  --print_angle_num_points=180
Read variable:  --load_input_file=1
Input file header: MC simulation: time=9999 bboxX=5 bboxY=5 bboxZ=5 beta=1e+08 kT=1e-08 V=125 U=-460.446 columns={type,x,y,z,ix,iy,iz,q0,q1,q2,q3,u}
atoms_to_load=160
Read variable:  --read_time_from_input=0
Note: Do not use time of input file header. Set read_time_from_input=1 to set time.
Read variable:  --read_bbox_from_input=0
Note: Do not use boundary box of input file header. Set read_bbox_from_input=1 to avoid this.
Read variable:  --read_kT_from_input=0
Note: Do not use kT from input file header. Set read_kT_from_input=1 to set kT to that of input header.
Note: Scale input coordinates (including volume) with coordinate_scaling, that is x_new = coordinate_scaling * x_old (coordinate_scaling=1.0  ->  no scaling).
Read variable:  --coordinate_scaling=1 
Read variable:  --type_of_volume_move=0
Read variable:  --delta_L_max=0.01 
Read variable:  --pressure=0 
Read variable:  --num_volume_moves_per_step=1
Note:  Set type_of_volume_move, delta_L_max, pressure to enable NpT simulations:
Note:      type_of_volume_move = 0   ->   No volume moves
Note:      type_of_volume_move = 1   ->   Isotropic volume moves                (cubic box       ; a*a*a ) 
Note:      type_of_volume_move = 2   ->   Anisotropic-orthorhombic volume moves (orthorhombic box; a*b*c ) 
Note:      type_of_volume_move = 3   ->   Anisotropic-tetragonal volume moves   (tetragonal box  ; a*a*c ) 
Note:      type_of_volume_move = 4   ->   Z-only moves                          (orthorhombic box; a*b*c ) 
Note: No barostat ( type_of_volume_move = 0 )
[Particle umbrella] = kappa*(N-N_{center})^2
umb_num_par sets the number of attempts pr MC step. Set to zero for no particle swat attempts.
Read variable:  --umb_num_par=0
Read variable:  --pin_particles=0
Note: Do not pin particles in space.
Read variable:  --max_step_size_translation=0.1 
Read variable:  --max_step_size_rotation=0.1 
Read variable:  --num_steps=100000
Read variable:  --moves_per_mc_step=1
Read variable:  --print_ener_freq=10
Read variable:  --print_conf_freq=100
Read variable:  --print_vmd_freq=10000
Note:   type_of_vmd_color=0   ->   particle- & arm energy
Note:   type_of_vmd_color=1   ->   orientational order
Read variable:  --type_of_vmd_color=0

   ..:: Simulation variables ::..
time = 0
beta = 1
boundary_box = [ 6 , 6 , 6 ]
seed = 1
skin_distance = 0.3
kT: 1
Potential energy: 26376.3
Volume: 216
Number of particles (sim): 160
Number of particles (potential): 160
Density: 0.740741
Average length of pair list: 11.7625

Begin run at 2024-10-03 09:39:05
Stat: step 0 (0%), t=0 U=26376.3 V=216 N=160 Q=0 UQ=0 Upin=0 moves=0 acc_trans=-nan%, acc_rot=-nan%, acc_vol=-nan%, acc_par=-nan%, nl_upd=1 t_wall=0.000397 t_est=inf vmd_prints=0 
Stat: step 2 (0.002%), t=2 U=13246.2 V=216 N=160 Q=0 UQ=0 Upin=0 moves=320 acc_trans=62.5%, acc_rot=100%, acc_vol=0%, acc_par=-nan%, nl_upd=3 t_wall=0.009363 t_est=468.2 vmd_prints=1 
Stat: step 4 (0.004%), t=4 U=8483.79 V=216 N=160 Q=0 UQ=0 Upin=0 moves=640 acc_trans=63.2812%, acc_rot=100%, acc_vol=0%, acc_par=-nan%, nl_upd=5 t_wall=0.013286 t_est=332.175 vmd_prints=1 
Stat: step 8 (0.008%), t=8 U=3176.26 V=216 N=160 Q=0 UQ=0 Upin=0 moves=1280 acc_trans=59.6094%, acc_rot=100%, acc_vol=0%, acc_par=-nan%, nl_upd=10 t_wall=0.021228 t_est=265.35 vmd_prints=1 
Stat: step 16 (0.016%), t=16 U=862.335 V=216 N=160 Q=0 UQ=0 Upin=0 moves=2560 acc_trans=59.1016%, acc_rot=99.9219%, acc_vol=0%, acc_par=-nan%, nl_upd=16 t_wall=0.035439 t_est=221.494 vmd_prints=1 
Stat: step 32 (0.032%), t=32 U=322.521 V=216 N=160 Q=0 UQ=0 Upin=0 moves=5120 acc_trans=62.0703%, acc_rot=99.8242%, acc_vol=0%, acc_par=-nan%, nl_upd=31 t_wall=0.064888 t_est=202.778 vmd_prints=1 
Stat: step 64 (0.064%), t=64 U=163.732 V=216 N=160 Q=0 UQ=0 Upin=0 moves=10240 acc_trans=67.0215%, acc_rot=99.5605%, acc_vol=0%, acc_par=-nan%, nl_upd=70 t_wall=0.125951 t_est=196.8 vmd_prints=1 
Stat: step 128 (0.128%), t=128 U=155.78 V=216 N=160 Q=0 UQ=0 Upin=0 moves=20480 acc_trans=71.8115%, acc_rot=99.4238%, acc_vol=0%, acc_par=-nan%, nl_upd=150 t_wall=0.251561 t_est=196.533 vmd_prints=1 
Stat: step 256 (0.256%), t=256 U=140.81 V=216 N=160 Q=0 UQ=0 Upin=0 moves=40960 acc_trans=74.5874%, acc_rot=99.3726%, acc_vol=0%, acc_par=-nan%, nl_upd=314 t_wall=0.500044 t_est=195.33 vmd_prints=1 
Stat: step 512 (0.512%), t=512 U=142.413 V=216 N=160 Q=0 UQ=0 Upin=0 moves=81920 acc_trans=75.8752%, acc_rot=99.165%, acc_vol=0%, acc_par=-nan%, nl_upd=621 t_wall=0.98557 t_est=192.494 vmd_prints=1 
Stat: step 1024 (1.024%), t=1024 U=131.037 V=216 N=160 Q=0 UQ=0 Upin=0 moves=163840 acc_trans=76.582%, acc_rot=99.1235%, acc_vol=0%, acc_par=-nan%, nl_upd=1267 t_wall=1.97503 t_est=192.874 vmd_prints=1 
Stat: step 2048 (2.048%), t=2048 U=136.869 V=216 N=160 Q=0 UQ=0 Upin=0 moves=327680 acc_trans=76.9156%, acc_rot=99.1232%, acc_vol=0%, acc_par=-nan%, nl_upd=2544 t_wall=3.94574 t_est=192.663 vmd_prints=1 
Stat: step 4096 (4.096%), t=4096 U=131.704 V=216 N=160 Q=0 UQ=0 Upin=0 moves=655360 acc_trans=77.028%, acc_rot=99.0919%, acc_vol=0%, acc_par=-nan%, nl_upd=5112 t_wall=7.89747 t_est=192.809 vmd_prints=1 
Stat: step 8192 (8.192%), t=8192 U=140.921 V=216 N=160 Q=0 UQ=0 Upin=0 moves=1310720 acc_trans=77.1263%, acc_rot=99.1056%, acc_vol=0%, acc_par=-nan%, nl_upd=10372 t_wall=15.9019 t_est=194.115 vmd_prints=1 
Stat: step 16384 (16.384%), t=16384 U=143.755 V=216 N=160 Q=0 UQ=0 Upin=0 moves=2621440 acc_trans=77.0485%, acc_rot=99.0857%, acc_vol=0%, acc_par=-nan%, nl_upd=20653 t_wall=31.7796 t_est=193.967 vmd_prints=2 
Stat: step 32768 (32.768%), t=32768 U=139.012 V=216 N=160 Q=0 UQ=0 Upin=0 moves=5242880 acc_trans=77.1053%, acc_rot=99.0922%, acc_vol=0%, acc_par=-nan%, nl_upd=41497 t_wall=63.7201 t_est=194.458 vmd_prints=4 
Stat: step 65536 (65.536%), t=65536 U=140.762 V=216 N=160 Q=0 UQ=0 Upin=0 moves=10485760 acc_trans=77.1105%, acc_rot=99.0933%, acc_vol=0%, acc_par=-nan%, nl_upd=82804 t_wall=127.412 t_est=194.415 vmd_prints=7 
Done with 100000 MC steps,
   corresponding to dt=100000.
End run at 2024-10-03 09:42:19
Run time of 194.341 seconds  = 3.23902 minutes  = 0.0539837 hours  = 0.00224932 days  = 0.000321332 weeks  = 0.00177363 months  = 0.000147803 years,  = 1.50603e-14 earth life-times, 
  corresponding to 514.558 MC steps per second,
  or 82329.3 particle MC steps per second.
Accepted 12339503 out of 16000000 attempted translation moves (77.1219%)
Accepted 15854785 out of 16000000 attempted rotation moves (99.0924%)
Accepted 0 out of 100000 attempted volume moves (0%)
Accepted 0 out of 0 attempted particle insertions/removals (-nan%)
Made 126449 neighbour list updates 
  corresponding to 1.26449 per MC step.
  corresponding to 0.00790306 per MC particle step.

   ..:: Simulation variables ::..
time = 100000
beta = 1
boundary_box = [ 6 , 6 , 6 ]
seed = 1
skin_distance = 0.3
kT: 1
Potential energy: 159.02
Volume: 216
Number of particles (sim): 160
Number of particles (potential): 160
Density: 0.740741
Average length of pair list: 10.275

Write final configuration to end.xyz and end.vmd.
Variable string: ../../bin/dna_colloid   --verbose=9 --seed=1 --kT=1 --chem_pot=0 --time=0 --bboxX=6 --bboxY=6 --bboxZ=6 --pair_eps=1 --pair_sigma=1 --arm_eps=1 --arm_r_center=1 --arm_r_width=0.2 --max_angle=40 --pair_12_eps=0 --pair_6_eps=0 --pair_4_eps=0 --cut_distance=1.2 --only_bind_to_lowest_arm_energy=0 --num_extra_sticky_spots=0 --type_of_orientational_order=1 --enable_orientational_order_threshold=0 --orientational_order_threshold_max=1.25 --orientational_order_threshold_min=0.75 --umbrella_orientational_order_kappa=0 --umbrella_orientational_order_center=0 --skin_distance=0.3 --num_added_particles=0 --print_random_numbers=0 --print_pair_num_points=250 --print_angle_num_points=180 --load_input_file=1 --read_time_from_input=0 --read_bbox_from_input=0 --read_kT_from_input=0 --coordinate_scaling=1 --type_of_volume_move=0 --delta_L_max=0.01 --pressure=0 --num_volume_moves_per_step=1 --umb_num_par=0 --pin_particles=0 --max_step_size_translation=0.1 --max_step_size_rotation=0.1 --num_steps=100000 --moves_per_mc_step=1 --print_ener_freq=10 --print_conf_freq=100 --print_vmd_freq=10000 --type_of_vmd_color=0
Happy ending of  ../../bin/dna_colloid
```

