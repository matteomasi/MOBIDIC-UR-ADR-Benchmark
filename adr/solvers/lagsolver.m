classdef lagsolver < handle
%% Semi-lagrangian solver for advection
% LAGSOLVER Implements semi-lagrangian integration schemes
% Dirichlet-type boundary condition at the inlet of the domain (left
% boundary x = 0)
% Neumann-type boudary condition at the outlet (right boundary)
%
%
% Usage. Create an instance:
%
%   obj = LAGSOLVER(x,C,v,dt)
%
%   obj = LAGSOLVER(x,C,v,dt,'solver','cubic_spline_advection') select
%       specific solver type. Default: cubic_spline_advection
%
%
%     x        Spatial coordinates (must be equally spaced)
%     v        Velocity 
%     D        Diffusion coefficient
%     dt       Time-step
%
%
% Available methods:
%   
%   LAGSOLVER.linear_advection(nt, C_bound) Linear interpolation
%   LAGSOLVER.cubic_advection(nt, C_bound) Cubic interpolation
%   LAGSOLVER.cubic_spline_advection(nt, C_bound) Cubic spline interpolation
%
% Preferred method: LAGSOLVER.cubic_spline_advection()
%
%
% Last update: 28/02/2024


properties
    x       % Spatial coordinate
    C       % Concentration 
    v       % Velocity
    dt      % Time-step
    dx      % Grid size
    solver  % Solver selection (Default: 'cubic_spline_advection')
end

methods
    %% Constructor
    function self = lagsolver(x,C_init,v,dt,varargin)
        self.x = x;
        self.C = C_init;
        self.v = v;
        self.dt = dt;
        self.dx = x(2)-x(1);
        
        % Parse optional input arguments
        p = inputParser;

        solverDefault = 'cubic_spline_advection';
        addParameter(p,'solver',solverDefault,@(x) ischar(x));

        parse(p,varargin{:});
        self.solver = p.Results.solver;
    end
 
    %% Method for setting new dt
    function setdt(self,dt)
        self.dt = dt;
    end

    %% Linear interpolation solver
    function C = linear_advection(self, nt, C_bound)
        % Propagate the current variable using a linear interpolation
        % nt : Number of time steps to be done
        % loop through time steps
        for t = 1:nt
            % temporary array
            y_temp = zeros(1,length(self.C));
            % loop through variable array
            for i = 1:length(self.C)
                % idx left to the departure point
                j = floor((self.x(i)-self.v*self.dt)/self.dx+1); % +1 because indexes start from 1
                % idx right to the departure point
                k = j+1;
                if j > 0
                    % linear interpolation
                    alpha = (self.x(i)-self.v*self.dt - self.x(j))/self.dx;
                    y_temp(i) = (1-alpha)*self.C(j) + alpha*self.C(k);
                else 
                    y_temp(i) = C_bound; % Left boundary
                end
            end
            % copy array to current time
            self.C = y_temp;
        end
        % Output the final state
        C = self.C;
    end
    
    
    %% Cubic interpolation solver
    function C = cubic_advection(self, nt, C_bound)
        % Propagate the current variable using a cubic interpolation scheme
        % nt : Number of time steps to be done  
        %loop through time steps
        for t = 1:nt
            % temporary array
            y_temp = zeros(1,length(self.C));
            % loop through array
            for i = 1:length(self.C)
                % idx left to departure point
                x_dep = self.x(i)-self.v*self.dt;
                j = floor(x_dep/self.dx)+1;
                if j > 1
                    % alpha
                    a = (self.x(i)-self.v*self.dt - self.x(j))/self.dx;
                    % calculate next time step
                    y_temp(i) = - a * (1-a)*(2-a)/6 * self.C(f(j-1));
                    y_temp(i) = y_temp(i) + (1-a^2)*(2-a)/2 * self.C(f(j));
                    y_temp(i) = y_temp(i) + a*(1+a)*(2-a)/2 * self.C(f(j+1));
                    y_temp(i) = y_temp(i) - a*(1-a^2)/6 * self.C(f(j+2));
                else
                    y_temp(i) = C_bound;
                end
            end
            self.C = y_temp;
        end
        function out = f(j)
           if j > length(self.C)
               out = rem(j,length(self.C));
           else
               out = j;
           end
        end
        % Output the final state
        C = self.C;
    end

    
    %% Cubic spline solver
    function C = cubic_spline_advection(self, nt, C_bound)
       
		% Propagate the current variable using a cubic spline interpolation
		% nt : Number of time steps to be done
		% loop through time steps
		% Code was vectorized for best performance
		for t = 1:nt
			% Cubic spline interpolation of current C values
			cs = spline(self.x, self.C); % Cubic spline coefficients
			shift = self.v*self.dt; % Spatial shift
			xi = self.x-shift; 
			k0 = xi<=0;  % Find points outside the domain
			xi(k0) = 0; % Set those points to left boundary
			yi = ppval(cs,xi); % Interpolated values
			yi(k0) = C_bound;  % Replace with left boundary
			self.C = yi;
		end
		% Output the final state
		C = self.C;

    end

    %% Solver helper function
    % Select solver based on user input, or use default solver
    function C = solve(self, nt, C_bound)

        switch self.solver
            case 'linear_advection'
                C = self.linear_advection(nt,C_bound);
            case 'cubic_advection'
                C = self.cubic_advection(nt,C_bound);
            case 'cubic_spline_advection'
                C = self.cubic_spline_advection(nt,C_bound);
            otherwise
                C = self.cubic_spline_advection(nt,C_bound);
        end

    end


end


end