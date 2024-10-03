# Model of DNA functionalized nano-particles

This program simulates a model of DNA functionalized nano-particles.
The program is written in C++, 
and developed on a Linux system. 

## Build from source

### Get the source code

Clone the repository hosted at github with

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

### Output files

