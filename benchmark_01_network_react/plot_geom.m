function plot_geom(node,ret,labels)
%% Plots the network's geometry
% labels:
%   0: no labels
%   1: node labels
%   2: reach labels
%
% Last revision: 28/02/2024

% Get node x,y and code
x_node = [node.xx];
y_node = [node.yy];
code_node = [node.code];

% Plot nodes
figure;
plot(x_node,y_node,'.')
hold on
axis equal

if labels == 1
    text(x_node,y_node,num2cell(code_node),'FontSize',8);
    t_title = 'Nodes';
end

% Number of reaches
N = length(ret);

% Plot network

for i = 1:N
    if labels == 2
        text(mean(ret(i).xx(1:2)),mean(ret(i).yy(1:2)),num2str(ret(i).code),'FontSize',8);
        t_title = 'Links';
    end
    plot(ret(i).xx(1:2),ret(i).yy(1:2),'-')
    hold on 
end
axis equal

xlabel('X (m)')
ylabel('Y (m)')
title(['Network topology - ' t_title])


end