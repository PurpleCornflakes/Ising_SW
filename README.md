# Potts Model with Swendsen-Wang Algorithm

## Potts Model

## Swendsen-Wang Cluster Algorithm

## Much More efficient than Metropolis Algorithm

## Code
1. iterate SW to a equilibrium state, output spin-lattice to a file 
2. read in the spin-lattice and continue iteration, record *energy* and *magnetization* at each step, check whether it's in equilibrium.
3. if in equilibrium, iterate and record slices at each step, draw MutualInfo versus distance.
