%% Rectangular box

clc
clear all
close all

% Example data

x=[0 5 5 0];
y=[0 0 5 5];

NodesBox=[x' y'];

StiffenerType=[1 2 NaN NaN];
StiffenerCenterDist=[0.9 0.5 NaN NaN];
DistanceEdgeMin=[0.5 0.5 NaN NaN];

ThicknessBox=[16 16 16 16]*1e-3;

StiffenerGeo{1}=[0 0 ; 83 275 ; 83+135 275 ; 300 0]*1e-3;
StiffenerGeo{2}=[0 0 ; 0 120 ; 30 200]*1e-3;

StiffenerThickness=[6 10]*1e-3;

opt=struct();
opt.plot=true
opt.plotstiffener=true
opt.plottext=true
opt.DistanceEdgeMin=DistanceEdgeMin

[Nodes,Elements,Thickness]=ClosedSectionStiffened(NodesBox,ThicknessBox,StiffenerType,StiffenerCenterDist,StiffenerGeo,StiffenerThickness,opt);

return
%% Sulafjorden

clc
clear all
close all

% Example data

x=[3.75 3.75+6 3.75+6+6.75 16.5-1.75 16.5-15.75 0];
y=[0 0 1.25 1.25+0.83 1.575+0.925 0.925];

NodesBox=[x' y'];

StiffenerType=[2 2 3 1 3 2];
StiffenerCenterDist=[0.9 0.9 0.6 0.6 0.6 0.9];
DistanceEdgeMin=[0.6 0.6 0.8 0.6 0.4 0.6];

ThicknessBox=[16 16 16 16 16 16]*1e-3;

StiffenerGeo{1}=[0 0 ; 83 275 ; 83+135 275 ; 300 0]*1e-3;
StiffenerGeo{2}=[0 0 ; 130 225 ; 130+190 225 ; 450 0]*1e-3;
StiffenerGeo{3}=[0 0 ; 100 150 ; 100+100 150 ; 300 0]*1e-3;

StiffenerThickness=[6 6 6]*1e-3;

opt=struct();
opt.plot=true
opt.plotstiffener=false
opt.plottext=true
opt.DistanceEdgeMin=DistanceEdgeMin

[Nodes,Elements,Thickness]=ClosedSectionStiffened(NodesBox,ThicknessBox,StiffenerType,StiffenerCenterDist,StiffenerGeo,StiffenerThickness,opt);

CrossxGenerateThinWalledSection(Nodes,Elements,Thickness,'TestSulafjorden','TestSulafjorden',true,'unit','m');

[A,Cx,Cy,Ix,Iy,Ixy,P]=PolygonMoments(Nodes(:,2:3))


%% Hålogland

clc
clear all
close all

x=cumsum([0 8 5.3 -1.7 -5.78 -(7.6+1.82) -1.7]);
y=cumsum([0 0 1.5 1.327 5.78*(3/100) -(7.6+1.82)*(3/100) -1.217]);
x=x-4;

NodesBox=[x' y'];

StiffenerType=[2 2 3 1 1 3 2];
StiffenerCenterDist=[0.875 0.95 0.35 0.6 0.6 0.35 0.95];
DistanceEdgeMin=nan*ones(size(StiffenerCenterDist));

ThicknessBox=[8 8 12 14 14 12 8]*1e-3;

StiffenerGeo{1}=[0 0 ; 83 275 ; 83+135 275 ; 300 0]*1e-3;
StiffenerGeo{2}=[0 0 ; 130 225 ; 130+190 225 ; 450 0]*1e-3;
StiffenerGeo{3}=[0 0 ; 0 150 ]*1e-3;
StiffenerThickness=[6 8 10]*1e-3;


DistanceEdgeStart=[0.5 0.475 1 0.68 0.3 0.5 1.2332];
N_stiff=[9 5 3 9 15 3 5]


opt=struct();
opt.plot=true
opt.plotstiffener=false
opt.plottext=true
opt.DistanceEdgeMin=DistanceEdgeMin
opt.DistanceEdgeStart=DistanceEdgeStart
opt.N_stiff=N_stiff

[Nodes,Elements,Thickness]=ClosedSectionStiffened(NodesBox,ThicknessBox,StiffenerType,StiffenerCenterDist,StiffenerGeo,StiffenerThickness,opt);


CrossxGenerateThinWalledSection(Nodes,Elements,Thickness,'TestHalogaland','TestHalogaland',true,'unit','m');

[A,Cx,Cy,Ix,Iy,Ixy,P]=PolygonMoments(Nodes(:,2:3))

%%
