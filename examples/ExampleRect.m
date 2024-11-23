%% Rectangular box

clc
clear all
close all

% Example data

x=[0 5 5 0];
y=[0 0 5 5];

nodes_box=[x' y'];

stiff_type=[1 2 NaN NaN];
stiff_cc=[0.9 0.5 NaN NaN];
distedgemin=[0.5 0.5 NaN NaN];

t_box=[16 16 16 16]*1e-3;

stiff_geo{1}=[0 0 ; 83 275 ; 83+135 275 ; 300 0]*1e-3;
stiff_geo{2}=[0 0 ; 0 120 ; 30 200]*1e-3;

stiff_t=[6 10]*1e-3;

opt=struct();
opt.plot=true
opt.plotstiffener=true
opt.plottext=true
opt.distedgemin=distedgemin

[Nodes,Elements,Thickness]=closedsectionstiffened(nodes_box,t_box,stiff_type,stiff_cc,stiff_geo,stiff_t,opt);

% CrossxGenerateThinWalledSection(Nodes,Elements,Thickness,'RectangleUnsym','examples\RectangleUnsym',true,'unit','m');

% [~,~,M_out]=CrossxExportParameters('examples\RectangleUnsymParameters.txt');

% [yc,zc,A,Iy,Iz,Iyz]=SectionParameters(Nodes,Elements,Thickness)

