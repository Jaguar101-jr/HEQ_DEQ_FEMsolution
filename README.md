# HEQ_DEQ_FEM_solver.
#This code is developed by HUSSEIN ABDUELHALEIM HUSSEIN MUHAMMED (https://github.com/Jaguar101-jr) as a tutorial to learn how to use FEM methods for solving PDEs.
# B.Sc.H. University of Khartoum (Sudan) 2015, M.Sc. from UPC-Qingdao (China) 2022.
# A Finite-Element method-based *matlab function [HEQ_DEQ_FEMsolution] to solve the  1D Heat Equation (Diffusion Equation) linear PDE.
a 1D domain with zero Dirichlet boundary conditions and an initial condition defined as u0 = sin(pi * x). 
The analytical solution for this problem can be found using the method of separation of variables.
The analytical solution for the heat equation with zero Dirichlet boundary conditions and the given initial condition is:
u(x, t) = exp(-pi^2 * alpha*t) * sin(pi*x)
where alpha is the thermal diffusivity.
You can extend the code to do: Convergence Study; Stability Check; Comparison with Other Solvers.
