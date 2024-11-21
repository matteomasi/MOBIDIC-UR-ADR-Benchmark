classdef diffsolver < handle
%% Collection of solvers for diffusion
% DIFFSOLVER Implements numerical schemes for the solution of diffusion equation
% Dirichlet-type boundary condition at the inlet of the domain (left
% boundary x = 0)
% Neumann-type boudary condition at the outlet (right boundary)
%
% Usage. Create an instance:
%
%   obj = DIFFSOLVER(x,C,D,dt) basic usage with step size dt
%
%   obj = DIFFSOLVER(x,C,D,dt,'solver','saulyev_solver_alt') select
%       specific solver type
%
%     x        Spatial coordinates (must be equally spaced)
%     C        Concentration 
%     D        Diffusion coefficient
%     dt       Time-step
%
% Optional properties (only for FD method)
%     phi      FD phi
%     theta    FD theta
%
%
% Available methods:
%   
%   DIFFSOLVER.pdepe_solver(nt, C_bound) uses pdepe matlab buit-in function
%   DIFFSOLVER.saulyev_solver_alt(nt, C_bound) uses two-direction Saul'yev method - uses mex file for increased speed
%
% Preferred method: DIFFSOLVER.saulyev_solver_alt()
%
%
% Last update: 28/02/2024


properties
    x       % Spatial coordinates
    C       % Concentration 
    D       % Diffusion coefficient
    dt      % Time-step
    dx      % Grid size
    phi     % FD phi
    theta   % FD theta
    solver  % Solver selection (Default: 'saulyev_solver_alt')
    F       % Fourier number (calculated from D, dt, dx)
end

methods
    %% Constructor
    function self = diffsolver(x,C_init,D,dt,varargin)
        self.x = x;
        self.C = C_init;
        self.D = D;
        self.dt = dt;
        self.dx = x(2)-x(1);
        self.phi = 0.5;
        self.theta = 0.5;
        
        % Parse optional input arguments
        p = inputParser;

        solverDefault = 'saulyev_solver_alt';
        addParameter(p,'solver',solverDefault,@(x) ischar(x));

        parse(p,varargin{:});
        self.solver = p.Results.solver;
        
        % Calculate Fourier number
        self.setdt(dt);
    end
    
    %% Method for setting new dt and recalculate Fourier number
    function setdt(self,dt)
        self.dt = dt;
        self.F = self.D*dt/(self.dx^2); % Callback update Fourier number
    end

    %% Matlab built-in solver (pdepe)
    function C = pdepe_solver(self, nt, C_bound)
        % nt : Number of time steps to be done
        % C_bound: Concentration at left boundary (inlet)
        % loop through time steps
        for k = 1:nt
            t_vec = linspace(0,self.dt,3);  % Local time vector (RK solver requires at least 3 sub-points)
            m = 0;   % Symmetry (0: no symmetry 1D problem)
            eqn = @(x,t,C,dCdx) diffpde(x,t,C,dCdx,self.D);  % Function handle of mainpde to pass extra parameters (D) 
            ic_fn = @(xi) ic(xi,self.x,self.C);
            bc_fn = @(xL,uL,xR,uR,t) bc(xL,uL,xR,uR,t,C_bound);
            sol = pdepe(m,eqn,ic_fn,bc_fn,self.x,t_vec); % SOLVE
            self.C = sol(end,:); % Take solution of the last sub-step
        end
        % Output the final state
        C = self.C;

        %% Diffusion PDE definition
        function [c,f,s] = diffpde(x,t,C,dCdx,D)
            c = 1;
            f = D*dCdx;
            s = 0;
        end

        %% Initial condition function
        function C0 = ic(xi,x,C_init)
            f = find(x==xi,1,'first');
            C0 = C_init(f);
        end
        
        %% Boundary condition function
        function [pL,qL,pR,qR] = bc(xL,uL,xR,uR,t,C_bound)
            pL = uL - C_bound; % Dirichlet at left boundary
            qL = 0;
            pR = uR;
            qR = 1; % Outflow right boundary
        end
    end
       
 

    %% Saulyev solver (integration in alternating directions)
    % With auto step size selector function
    function C = saulyev_solver_alt(self, nt, C_bound)
        
        dt = self.dt;
		theta = self.D*dt/(self.dx^2);
				
		% Assign current C state as initial condition
		C_init = self.C;
		CLR = self.C; CRL = self.C;
		
		% Loop through time steps
		for k = 1:nt
			
			% A) L-R direction
			for i = 1:length(CLR)
				if i == 1 % left boundary
					solA = theta*C_bound;
				else
					solA = theta*CLR(i-1);
				end
				
				solB = (1-theta)*C_init(i);
				
				if i < length(CLR)
					solC = theta*C_init(i+1);
				else
					solC = theta*C_init(i);
				end

				% L-R Solution
				CLR(i) = (solA+solB+solC)/(1+theta);
			end
			
			% B) R-L direction
			for i = length(CRL):-1:1
				if i == length(CRL) % right boundary (take from LR solution)
					solA = theta*CLR(end);
				else
					solA = theta*CRL(i+1);
				end
				
				solB = (1-theta)*C_init(i);
				
				if i > 1
					solC = theta*C_init(i-1);
				else
					solC = theta*C_init(i);
				end
				
				% R-L Solution
				CRL(i) = (solA+solB+solC)/(1+theta);
 
			end
			
			% Average L-R and R-L solutions
			C = (CLR+CRL)/2;
			
			% Update initial condition for next step k
			C_init = C;
			
		end
		% Update to the final state
		self.C = C;
             
    end
    

    %% Solver helper function
    % Select solver based on user input, or use default solver
    function C = solve(self, nt, C_bound)

        switch self.solver
            case 'pdepe_solver'
                C = self.pdepe_solver(nt,C_bound);
            case 'saulyev_solver_alt'
                C = self.saulyev_solver_alt(nt,C_bound);
            otherwise
                C = self.saulyev_solver_alt(nt,C_bound);
        end

    end
    
end


end