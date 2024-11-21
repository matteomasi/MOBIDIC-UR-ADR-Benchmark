%% HFCW 1D - Horizontal Flow Constructed Wetland Model - Tracer Validation
% Simplified 1D advection/dispersion model 
% CWM1 - Constr. Wetland Mod. 1 - Langergraber et al. (2009)
%
% Tracer data from Boog et al. (2019)
% Wetland construction: Nivala et al. (2013)
%
% Wetland data:
%   4.7m in length x 1.2m in width with a saturated depth of 0.9 m.
%   Hydraulic loading rate (HLR): 0.576 m3/d
%   Flowrate Q = 6.667 E-6 m3/s = 0.4 L/min
%	12 L of primarily treated domestic sewage every 30 min at a dosing rate of 5 L/min
%
%   URANINE (Fluorescein)
%   Tracer concentration inlet: 2.5 mL of a 200 g/L (0.5g) solution in 12 L  ==>    Ctracer = 41.66 mg/L  in the 12L tank
%   Tracer mass recovery at the outlet: 92% (uranine, control test)
%
%   BROMIDE Br-
%   60 g of dried KBr -> 40 g of Br-
%   40 g Br- in 12 L -> 3.33 g/L
%
% Last revision: 28/02/2024
clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL SETUP

% Load parameters, boundary conditions, performs model initialization
parameters

% Add solvers/helpers path
addpath('../adr/solvers')
addpath('../adr/helpers')


% Initialize solver objects
for s = 1:Ncomp
    s_diff(s) = diffsolver(x,C_init(s,:),D,dt,'solver',diffSolver);
    s_adv(s)  = lagsolver(x,C_init(s,:),v,dt,'solver',advSolver);
end


%% MAIN LOOP
% Solve using coupling helper
a = tic;
[tout,C,stats] = coupling(s_diff,s_adv,C_L,'dtout',dt_output);
tot = toc(a);
disp(['Total diffusion computation (s): ' num2str(stats.timer(1))])
disp(['Total advection computation (s): ' num2str(stats.timer(2))])
disp(['Total computation time (s): ' num2str(tot)])


%% Display PLOTS (see plots.m for configuration)
plots


