clear all;
close all;
clc;
%% Ploting

load pos_nodes.dat
load vertiic.dat
load E_field.dat
load H_field.dat

N_node = size(pos_nodes,1);
N_patch = size(vertiic,1);

figure
patch('Faces',vertiic ,'Vertices',pos_nodes, 'FaceColor', 'b');
axis equal
xlabel('x(m)');
ylabel('y(m)');
hold on

scatter(E_field(:,1),E_field(:,2),10,'red')
hold on

scatter(H_field(:,1),H_field(:,2),10,'green')
hold off
