%% PLOTS
% BENCHMARK TEST - WASP 8
% Run this file after the end of the simulation (main.m)
%
% Last revision: 28/02/2024

%%%%%% PARAMETERS %%%%%%

% Node(s) to be plotted
node_plot = [2,4,6];

% Reach(es) to be plotted
reach_plot = [5];

% Number of time subdivisions (max)
Nt_sub = 35;

%%%%%%%%%%%%%%%%%%%%%%%%

%% Load WASP data
wasp = readtable('wasp8_output.xlsx');
t_wasp = minutes(wasp.E1_t-wasp.E1_t(1));
reach_wasp = {'C1','E1','E20'};
node_wasp = [2,4,6];

%% Plot network
plot_geom(node,ret,1);
plot_geom(node,ret,2);


%% Plot concentration at selected nodes
Csel = zeros(length(Cnodes),length(node_plot));
for i = 1:length(Cnodes)
    for j = 1:length(node_plot)
        k = find([node.code]==node_plot(j), 1, 'first');
        Csel(i,j) = Cnodes{i}(k);
    end
end


figure;
lgnd = {};
for i = 1:size(Csel,2)
    plot(params.t,Csel(:,i))
    hold all
    lgnd{end+1} = ['Node ' num2str(node_plot(i))];
end

for i = 1:length(reach_wasp)
    var_str = [reach_wasp{i} '_C'];
    plot(t_wasp, wasp.(var_str),'o')
    lgnd{end+1} = ['WASP Node ' num2str(node_wasp(i))];
end
    
legend(lgnd)
xlabel('Elasped time (min)')
ylabel('Concentration (g/m^3)')
xlim([0 200])

