%% HFCW 1D - Horizontal Flow Constructed Wetland Model - PLOTS
% Simplified 1D advection/dispersion model coupled to CWM1 model
% CWM1 - Constr. Wetland Mod. 1 - Langergraber et al. (2009)
%
%
% Last revision: 28/02/2024


%% Load tracer data
data_Br = csvread('Boog_et_al_2019_data/Boog_et_al_2019_Br.csv');
data_Ur = csvread('Boog_et_al_2019_data/Boog_et_al_2019_uranine.csv');
t_Ur_i = linspace(0,15,50);
Ur_i = interp1(data_Ur(:,1),data_Ur(:,2),t_Ur_i,'linear','extrap');
Ur_i = Ur_i - data_Ur(1,2);
% Boog model
model_Br = csvread('Boog_et_al_2019_data/Boog_et_al_2019_Br_model.csv');
model_Ur = csvread('Boog_et_al_2019_data/Boog_et_al_2019_uranine_model.csv');
t_mod = linspace(0,15,100);
Br_mod = interp1(model_Br(:,1),model_Br(:,2),t_mod,'linear','extrap');
Ur_mod = interp1(model_Ur(:,1),model_Ur(:,2),t_mod,'linear','extrap');

%% Plot tracer validation
t_out = (0:dt_output:t(end))'/3600/24;

% Br-
figure;
plot(data_Br(:,1),data_Br(:,2),'o','MarkerSize',9,'Color','black')
hold all
plot(t_mod,Br_mod,'--k','LineWidth',1.5)
plot(t_out,C(:,end,1),'LineWidth',1.5,'Color',[0.2 0.7 0.2])
xlabel('Elapsed time (days)')
ylabel('Concentration (mg/L)')
legend('Boog et al. 2019 - Br^- data','Boog et al. 2019 - Br^- model','MobidicU-ADR - Br^-')

% Uranine
figure;
plot(t_Ur_i,Ur_i,'o','MarkerSize',9,'Color','black')
hold all
plot(t_mod,Ur_mod,'--k','LineWidth',1.5)
plot(t_out,C(:,end,2),'LineWidth',1.5,'Color',[0.2 0.7 0.2])
xlabel('Elapsed time (days)')
ylabel('Concentration (mg/L)')
legend('Boog et al. 2019 - Fluorescein data','Boog et al. 2019 - Fluorescein model','MobidicU-ADR - Fluorescein')

%% Mass balance tracer
M_recover = sum(C(:,end,1))*dt_output*Q*1000;    % Recovered at outlet
M_inlet = sum(C_L(:,1))*dt*Q*1000;           % Calculated from boundary condition
disp(['Br- mass balance: ' num2str(M_recover/1000) ' g (inlet mass: ' num2str(M_inlet/1000) ' g). Mass balance error: ' num2str(100-M_recover/M_inlet*100) '%'])

M_recover2 = sum(C(:,end,2))*dt_output*Q*1000;    % Recovered at outlet
M_inlet2 = sum(C_L(:,2))*dt*Q*1000;           % Calculated from boundary condition
disp(['Fluorescein mass balance: ' num2str(M_recover2) ' mg (inlet mass: ' num2str(M_inlet2) ' mg). Mass balance error: ' num2str(100-M_recover2/M_inlet2*100) '%'])




