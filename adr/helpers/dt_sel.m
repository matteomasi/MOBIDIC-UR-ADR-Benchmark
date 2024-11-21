function [dtn, dx_adj] = dt_sel(dt, dx, v, DL)
%% Selection of integration step
% Automated selection of integration time steps
%
%
% INPUT:
%   dt: main local integration time-step in seconds 
%   dx: initial dx assignments (dx_max), in meters
%   v: velocity (m/s)
%   DL: dispersion coefficient (m2/s)
%
% OUTPUT
%   dtn: number of substeps 
%   dx_adj: adjusted dx
%
% Last update: 28/02/2024

% Assign dtn

Cr_max = 5;
dt_opt = Cr_max.*dx./v;
dtn = ceil(dt./dt_opt);


% Assign dx 
dx_adj = ones(size(v)).*dx;


end