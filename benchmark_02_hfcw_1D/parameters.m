%% HFCW 1D - PARAMETERS and BOUNDARY CONDITIONS
% Simplified 1D advection/dispersion model 
%
% Last revision: 05/05/2022
% M.M.2022

%% PARAMETERS
%%%%%%%%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_fin = 25;       % Simulation duration (days)
dt = 100;          % Transport time step (s)

dt_output = 1000; % Output time step (s)


%%%%%%%%%%%%%%% TRANSPORT MODEL PARAMETERS %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 4.7;        % Domain length (m)
nx = 50;        % Number of grid points
D = 1.6e-5;     % Dispersion coefficient (m2/s).  D = alpha * v = 1.54e-5 m2/s   (alpha = 2.5, dispersion length fast flow domain)
Q = 6.667e-6;   % Flowrate (m3/s)
w = 1.2;        % Width (m)
d = 0.9;        % Depth (saturated) (m)

diffSolver = 'saulyev_solver_alt';      % Diffusion solver
advSolver = 'cubic_spline_advection';   % Advection solver


%%%% BOUNDARY CONDITIONS
% Step boundary condition
tb_0 = 0;               % Step concentration start (s)
tb_1 = 1800;            % Step concentration stop (s) 30 min dosing event
C_Ur_L = 500/12;         % Boundary concentration uranine (mg/L)
C_Br_L = 40/12*1000;    % Boundary concentration Br- (mg/L)

%%%% TRANSPORT MODEL
x = linspace(0,L,nx); % Array of grid points (MUST be equidistant!)
dx = x(2) - x(1); % spatial grid point distance 
v = Q/w/d;      % Advective velocity (m/s)
% Courant, Peclet and Fourier numbers
Cr = v*dt/dx;
Pe = v*dx/D;
F = D*dt/(dx^2);
disp(['Courant = ' num2str(Cr)])
disp(['Peclet = ' num2str(Pe)])
disp(['Fourier = ' num2str(F)])

%%%% MODEL INITIALIZATION
t_fin_s = t_fin*24*3600;            % Simulation duration (s)
t = 0:dt:t_fin_s;                   % Time vector (s)
nt = length(t);                     % Number of time steps

Ncomp = 2;                          % Number of model components 
C_init = zeros(Ncomp,nx);           % Initial condition (zero)
C_L = zeros(nt, Ncomp);             % Zero boundary condition
C_L(t>tb_0&t<=tb_1,1) = C_Br_L;     % Step boundary condition Br-
C_L(t>tb_0&t<=tb_1,2) = C_Ur_L;     % Step boundary condition uranine


