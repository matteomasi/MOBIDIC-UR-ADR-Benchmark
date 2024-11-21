%% Network advection-dispersion solver
% BENCHMARK TEST - WASP 8
%
% Last revision: 28/02/2024


clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.file.hydraulics = 'network';         % Input .mat file containing Qret, ret, node 
params.file.C_input    = 'input';           % Input .mat file containing Cin


% Dispersion coefficient
params.D = 1;     % m2/s

% First order reaction coefficient
params.K = 10;    % 1/days


% Solver configuration
params.solver.type = 'Lagrangian';
params.solver.lag_adv = 'cubic_spline_advection';
params.solver.lag_diff = 'saulyev_solver_alt';


% Global time-stepping settings (hydraulics)
params.dt = 0.5;         % Time step (minutes) 0.5
params.tfin = 300;        % Simulation duration (minutes)  80


% Local time-stepping
params.dt_local = 0.1;      % Initial local time-step (minutes) Lag: 0.1, FD: 0.005

% Local numerical grid
params.dx_max = 1;        % Maximum distance between points (m) 0.25

% Initial concentration in channels
params.Cinit = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
% Temporarily add the directory with solvers and helpers to Matlab path 
addpath('../adr/solvers')
addpath('../adr/helpers')

% Load the network/hydraulics file
load([params.file.hydraulics '.mat']);

% Load the pollutant input file
load([params.file.C_input '.mat']);

% Number of nodes/reaches
Nnode = length(node);
Nreach = length(Qret);

% Number of global time steps
params.Nt = floor(params.tfin/params.dt)+1;
params.t = 0:params.dt:params.tfin;

% Number of local time steps
params.Nt_local = ceil(params.dt/params.dt_local);    


%% Dispersion coefficients
DL = params.D*ones(Nreach,1);
tmpcell = num2cell(DL);
[Crettmp(1:Nreach).DL] =  tmpcell{:};
% Assign velocity 
tmpcell = num2cell([Qret.Vaverage]);
[Crettmp(1:Nreach).v] = tmpcell{:};
% Copy to all time steps
Cret = cell(params.Nt,1);
for i = 1:params.Nt
    Cret{i} = Crettmp;
end

clear Crettmp


%% Calculate inflow/outflow to/from nodes

Qin = zeros(Nnode,1); Qout = Qin;
for i = 1:Nnode
    % Calculate node inflow flow rate Qin
    Qin(i) = 0;
    if ~isempty(node(i).linkIN)
        Qin(i) = sum( [Qret(node(i).linkIN).Qout] );
    end
    
    % Calculate node outflow rate Qout
    Qout(i) = 0;
    if ~isempty(node(i).linkOUT)
        Qout(i) = sum( [Qret(node(i).linkOUT).Qout] );
    end
    
    % Exception: OUTLET
    if strcmpi(node(i).type, 'outlet')
        Qout(i) = sum( [Qret(node(i).linkIN).Qout] );
    end
end



%% Calculate initial concentrations Cnodes{1} in the nodes from Cbound{1}
% Initialize variable that stores the concentration at the nodes at each t
Cnodes = cell(params.Nt,1);
Cnodes{1} = zeros(Nnode,1);
for i = 1:Nnode
    if ~isnan(Cbound{1}(i).Cnode) 
        Cnodes{1}(i) = Cbound{1}(i).Cnode;  % Fixed concentration
    else
        Cnodes{1}(i) = Cbound{1}(i).Fnode/Qout(i); % Calc. from influx
    end
end


%% Selection of the integration parameters 
% Vector of channel water velocities
v = [Qret.Vaverage];

% Local time-step, conversion to seconds
params.dt_local_s = params.dt_local*60; 

% Selection of local integration time steps
[params.dtn_local, params.dx_adj] = dt_sel(params.dt_local_s,params.dx_max,v,DL);
params.dt_local_adj = params.dt_local_s./params.dtn_local;

% Calculate characteristic numbers
adim = calc_numbers(params.dt_local_adj, params.dx_adj, v, DL);
params.Cr = adim.Cr; 
params.Pe = adim.Pe; 
params.F = adim.F;



% Vector of channel lengths
L = [ret.L];



% Assign numerical model properties to each channel
Cprop = cell(Nreach,1);
for i = 1:Nreach
    % Local space vectors (selection of dx_adj for each single channel still to
    % be implemented in dt_sel() function)
    Cprop{i}.num_x = ceil(L(i)/params.dx_adj(i)); 
    Cprop{i}.x = linspace(0,L(i),Cprop{i}.num_x);
    Cprop{i}.dx = Cprop{i}.x(2)-Cprop{i}.x(1);

    % Local time vectors (substeps)
    Cprop{i}.dt = params.dt_local_adj(i);
    Cprop{i}.num_t = params.dtn_local(i);
    Cprop{i}.t = linspace(0,Cprop{i}.dt,Cprop{i}.num_t);

    % Calculate Cr, Pe and F
    adim = calc_numbers(Cprop{i}.dt, Cprop{i}.dx, v(i),DL(i));
    Cprop{i}.Cr = adim.Cr;Cprop{i}.Pe = adim.Pe;Cprop{i}.F = adim.F;

    % Misc
    Cprop{i}.L = L(i);
    Cprop{i}.v = v(i);
    Cprop{i}.D = DL(i);


    % Initialize concentrations (global time step)
    Cret{1}(i).C = zeros(size(Cprop{i}.x));
 
end
    

% Initialize local Cret
Cret_local = cell(params.Nt_local+1,1);

% Initialize local Cnodes
Cnodes_local = cell(params.Nt_local+1,1);


%% MAIN LOOP
ticAll = tic;
fprintf('-------------------------------------------------\n')
fprintf('Simulation started %s \n',datestr(now))
fprintf('-------------------------------------------------\n')

for k = 2:params.Nt % Global time step
    ticPartial = tic;
    fprintf('Time step index %4.0f, Simulation time: %0.4f min...',k,params.t(k));
    
    Cret_local{1} = Cret{k-1};      % Initial conditions
    Cnodes_local{1} = Cnodes{k-1};  % Boundary conditions
    
    % Iteration over local time steps
    for kL = 2:params.Nt_local+1

        % Update node concentrations
        Cnodes_local{kL} = zeros(Nnode,1);
        for i = 1:Nnode
            if ~isnan(Cbound{k-1}(i).Cnode) 
                Cnodes_local{kL}(i) = Cbound{k-1}(i).Cnode;  % Fixed concentration
            else
                % Upstream reaches inflow
                if ~isempty(node(i).linkIN)
                    Cend = [];
                    for ii = 1:length(node(i).linkIN)
                        Cend = [Cend Cret_local{kL-1}(node(i).linkIN(ii)).C(end)]; % Take the concentration at the reach's outlet
                    end
                    Fu =  sum( [Qret(node(i).linkIN).Qout] .* Cend );
                else
                    Fu = 0;
                end

                Cnodes_local{kL}(i) = (Fu + Cbound{k-1}(i).Fnode)/Qout(i); % Calc. from influx
            end
        end
        
        % Iteration over all the reaches
        Cret_kLm1 = Cret_local{kL-1}; % This is necessary for parfor
        Cret_kL = Cret_local{kL-1};   % Initialize variable
        
        for i = 1:Nreach
            C_L = Cnodes_local{kL}(ret(i).n1); % Boundary condition upstream
            C_init = Cret_kLm1(i).C;       % Initial concentration from the previous step
            v = Cprop{i}.v; D = Cprop{i}.D; x = Cprop{i}.x; dt = Cprop{i}.dt; nt = Cprop{i}.num_t;

            % Initialize solver(s)
            if strcmp(params.solver.type, 'Lagrangian') 
                % LAGRANGIAN
                s1 = lagsolver(x,C_init,v,dt);
                s2 = diffsolver(x,C_init,D,dt);
                for jj = 1:nt
                    s2.C = s1.(params.solver.lag_adv)(1,C_L);   % Advection
                    s1.C = s2.(params.solver.lag_diff)(1,C_L);  % Diffusion
                    s1.C = s1.C.*exp(-params.K*dt/3600/24);     % First-order reaction
                    
                end
            else
                % FULL
                s1 = fullsolver(x,C_init,v,D,dt);
                s1.C = s1.(params.solver.full)(nt,C_L,params.solver.fd_scheme,0.5); % FD
                s1.C = s1.C.*exp(-params.K*dt/3600/24);             % First-order reaction
            end

            % Update state for i-th reach
            Cret_kL(i).C = s1.C;   
        end
        
        % Save current state
        Cret_local{kL} = Cret_kL;        
        
    end
    
    % Update local concentrations and node concentrations
    for h = 1:Nreach
        Cret{k}(h).C = Cret_local{end}(h).C;
    end
    Cnodes{k} = Cnodes_local{end};
    
   
    
    fprintf('Elasped time %0.2f s. ',toc(ticPartial));
    fprintf('Percent completed: %.1f %% \n',100*k/(params.Nt));

end



fprintf('--------\n')
timeTot = toc(ticAll);
fprintf('Computation time\t\t\t %s\n ', datestr(seconds(timeTot),'HH:MM:SS') )
fprintf('-------------------------------------------------\n')


%% Calculate mass balance error
Fin_all = zeros(params.Nt,Nnode); % g/s
Min = zeros(Nnode,1);
for k = 1:params.Nt
    for i = 1:Nnode
        if ~isnan(Cbound{k}(i).Cnode) 
            Fin_all(k,i) = Cbound{k}(i).Cnode*Qout(i); % Fixed concentration
        else
            Fin_all(k,i) = Cbound{k}(i).Fnode; % Calc. from influx
        end
    end
end
for i = 1:Nnode
    % Find outlet node
    if strcmpi(node(i).type, 'outlet')
        outlet_node = i;
    end
end

Coutlet = zeros(length(Cnodes),1);
for i = 1:length(Cnodes)
    Coutlet(i) = Cnodes{i}(outlet_node);
end


%% Display PLOTS (see plots.m for configuration)
plots


