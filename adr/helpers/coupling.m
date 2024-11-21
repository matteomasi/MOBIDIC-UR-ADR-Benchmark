function [tout,C,diagnostics] = coupling(diff,adv,C_L,varargin)
% SNIA coupling scheme (SNIA: sequential non iterative approach)
% The initial condition is taken from the advection object C0 = adv.C
% The initial step size is taken from adv.dt
% The simulation duration is calculated from the length of C_L
%
% Objects (diff, adv) are passed by reference (they are not duplicated
% within this function)
%
% Usage:
%   COUPLING(diff,adv,C_L)                  Basic usage only transport coupling
% 
% Options:
%   'dtout' (float). Resample output vectors with dtout interval
%
%
% Inputs:
%   diff         Single object or vector of objects diffusion solver class
%   adv          Single object or vector of objects advection solver class
%   C_L          Time-dependent inlet boundary condition (nt*Ncomp)
%
% Output:
%   t_out               Time vector
%   C                   Concentration vector (nt*Ncomp)
%   diagnostics.timer   Execution elapsed time [diff adv react] (s)
% 
%
% References Richardson extrapolation (usually for adaptive time stepping):
%   Bell and Binning (2004)
%   Fahs et al. (2008) 
%
%
% Last update: 28/02/2024


%% PARAMETERS


%% Input parsing
% Check object input size
n_diff = numel(diff);
n_adv  = numel(adv);
if n_diff ~= n_adv
    error('Error: number of elements in input vectors are not equal')
end
% Check dt sizes
dt1 = zeros(n_diff,1);
dt2 = zeros(n_adv,1);
for i = 1:n_diff
    dt1(i) = diff(i).dt;
    dt2(i) = adv(i).dt;
end

if ~all(dt1==dt2)
    error('time steps (dt) of input objects are not equal')
end

% Store user dt
dt = dt1(1);

% Parse varargin
p = inputParser;

dtoutDefault = dt;
addParameter(p,'dtout',dtoutDefault,@(x) isnumeric(x));

parse(p,varargin{:});
dtout = p.Results.dtout;

%% MAIN LOOP
Ncomp = n_adv;
nx = length(adv(1).x);      % Number of mesh elements
nt = size(C_L,1);           % Number of steps with fixed size dt
t = 0:dt:(nt-1)*dt;         % Time vector
tout = (0:dtout:t(end))';   % Downsampled time vector
ntout = length(tout);       % Number of elements tout
C = zeros(ntout,nx,Ncomp);  % Create array for storing evolution of C(t,x,comp)
diagnostics.timer = zeros(1,3);         % Execution time [adv diff react]


% Misc variables initialization
kout = 2;
a_tm = 0; d_tm = 0; r_tm = 0; % Initialize timers

for k = 2:nt
    for i = 1:Ncomp
        tic; adv(i).C = diff(i).solve(1,C_L(k,i)); d_tm = toc+d_tm;
        tic; diff(i).C = adv(i).solve(1,C_L(k,i)); a_tm = toc+a_tm;
        
        if rem(t(k),dtout) == 0
            % Concentration
            C(kout,:,i) = diff(i).C;
            if i == Ncomp
                kout = kout+1;
            end
        end
        
    end
end


% Output diagnostics
diagnostics.timer = [d_tm a_tm r_tm];
    
end