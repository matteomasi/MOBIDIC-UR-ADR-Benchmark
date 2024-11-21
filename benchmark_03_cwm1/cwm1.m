function C = cwm1(t_fin, params, init_cond)
%% CWM1 handling function
% Build the ODE system and run the ODE solver 
% Input:
%   - t_fin
%   - params: vector of model paramenter
%   - init_cond: matrix [NxM], where N is the number of instances
%     (transport model mesh), and M is the number of components
%
% Output:
%   - C: Concentration matrix [NxM] 
%
% Last update: 28/02/2024

    ode_eqn = @(t,C) cwm1_odesystem(t,C,params); % Function handle system of equation

    C = zeros(size(init_cond));
    for i = 1:size(init_cond,1)
        [~,Ctmp] = ode23(ode_eqn, [0,t_fin] , init_cond(i,:)); % Solve
        C(i,:) = Ctmp(end,:);
    end
    
end