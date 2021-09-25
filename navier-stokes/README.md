# Navier-Stokes Method

## Compilation

```
gfortran -o navier-stokes_bin navier-stokes.f90  -O2
```

## Configure the input file

**input.txt** - contains parameters of the task
	
EXAMPLE

```
0.0	0.1	101	! START_X, END_X, NUMB_X
0.0	0.05	101	! START_Y, END_Y, NUMB_Y
0.01	0.0		! U0, V0
0.01			! P0
1.E-6			! MU
0.5			! CFL
100000			! MAX_ITER
```

*NOTE:* symbols after the sign "!" are ignored

## Launch

```
./navier-stokes_bin
```

## Output

**output.plt** - 2D - boundary layer field
**residuals.dat** - residuals
**cf.dat** - skin friction coefficient on the bottom wall