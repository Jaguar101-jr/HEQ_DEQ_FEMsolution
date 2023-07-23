
%% This code is developed by HUSSEIN ABDUELHALEIM HUSSEIN MUHAMMED (https://github.com/Jaguar101-jr)
%% As a tutorial to learn how to use FEM methods for solving PDEs.
%% B.Sc.H. University of Khartoum (Sudan) 2015, M.Sc. from UPC-Qingdao (China) 2022.
%% A Finite-Element method-based *matlab function [HEQ_DEQ_FEMsolution] to solve the  1D Heat Equation (Diffusion Equation) linear PDE.
%% a 1D domain with zero Dirichlet boundary conditions and an initial condition defined as u0 = sin(pi * x). 
%% The analytical solution for this problem can be found using the method of separation of variables.
%% The analytical solution for the heat equation with zero Dirichlet boundary conditions and the given initial condition is:
%% u(x, t) = exp(-pi^2 * alpha*t) * sin(pi*x)
%% where alpha is the thermal diffusivity.
%% You can extend the code to do: Convergence Study; Stability Check; Comparison with Other Solvers.


%% Define the function
function FEM_Heat_Equation_1D()
    % Parameters
    L = 1;          % Length of the domain
    nx = 50;        % Number of elements (reduced for demonstration)
    dx = L / nx;    % Element size
    nt = 100;       % Number of time steps (reduced for demonstration)
    dt = 0.001;     % Time step size
    alpha = 0.01;   % Thermal diffusivity

    % Initial condition
    u0 = @(x) sin(pi * x);

    % Create the mesh
    x = linspace(0, L, nx + 1); % Node coordinates
    elements = [(1:nx)', (2:nx+1)']; % Element connectivity

    % Assemble the stiffness matrix and load vector
    K = zeros(nx+1, nx+1);
    F = zeros(nx+1, 1);
    %% iterative FEM solution
    for i = 1:nx
        % Element matrix
        ke = element_matrix(x(i), x(i+1), alpha);

        % Assemble into the global stiffness matrix
        K(elements(i, :), elements(i, :)) = K(elements(i, :), elements(i, :)) + ke;

        % Element load vector
        fe = element_load(x(i), x(i+1), u0);

        % Assemble into the global load vector
        F(elements(i, :)) = F(elements(i, :)) + fe;
    end

    % Apply zero Dirichlet boundary conditions
    K(1, :) = 0; K(1, 1) = 1;
    K(nx+1, :) = 0; K(nx+1, nx+1) = 1;
    F(1) = 0;
    F(nx+1) = 0;

    % Time-stepping loop
    u = u0(x); % Initial condition
    for t = 1:nt
        u = u + dt * (K \ F);
    end

    % Analytical solution
    t_final = nt * dt;
    u_analytical = exp(-pi^2 * alpha * t_final) * sin(pi * x);

    % Plot the numerical and the analytical solutions
    figure;
    plot(x, u, '-o', 'DisplayName', 'Numerical Solution');
    hold on;
    plot(x, u_analytical, 'r-', 'DisplayName', 'Analytical Solution');
    xlabel('Displacement (cm) X-axis');
    ylabel('Temperature (degree)');
    title('1D Heat Equation FEM Solution');
    legend('Location', 'northwest', 'FontSize', 7, 'NumColumns', 2);
    hold off;
end

function ke = element_matrix(x1, x2, alpha)
    % Element stiffness matrix for 1D linear element
    ke = alpha / (x2 - x1) * [1, -1; -1, 1];
end

function fe = element_load(x1, x2, u0)
    % Element load vector for 1D linear element
    fe = (x2 - x1) / 2 * [u0(x1) + u0(x2)];
end
