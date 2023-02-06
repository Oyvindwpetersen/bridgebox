%%

% clc
% clear all
% close all
% 
% Folder=[BaseFolder '\' 'Abaqus\Suspensionbridge\crosssection_prop\'];
% t=20e-3
% E=210e9
% v=0.3
% 
% 
% Nodes=[0 0 ; 2 0 ; 2 1 ; 1 3 ; 0 1];
% 
% % x=[7.5 7.5+12 7.5+12+13.5 7.5+12+13.5-3.5 1.5 0]
% % y=[0 0 2.5 4 6 2]
% % Nodes=[x' y'];
% 
% [yc,zc,ys,zs,A,Iyy,Izz,Iyz,It,theta,P]=abdbeam_csprop(t,E,v,Nodes,Folder,'offsetgeometry',false);


%%

clc
clear all
close all

Folder=[BaseFolder '\' 'Abaqus\Suspensionbridge\crosssection_prop\'];
E=210e9
v=0.3

% Nodes=[0 0 ; 2 0 ; 2 1 ; 1 3 ; 0 1];
% Elements=[];
% t_mat=[ [1:5].' [0.02 0.01 0.02 0.02 0.01].' ];

Nodes=[0 0 ; 2 0 ; 2 1 ; 1 3 ; 0 1];
Nodes=[ [1:5]' Nodes];
Elements=[1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 1 ; 3 5];
Elements=[ [1:6]' Elements];
t_mat=[ [1:6].' [0.02 0.01 0.02 0.02 0.01 0.03].' ];

[yc,zc,ys,zs,A,Iyy,Izz,Iyz,It,theta,P]=abdbeam_csprop(E,v,Nodes,Elements,t_mat,Folder,'offsetgeometry',false);
