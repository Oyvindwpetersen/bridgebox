%%

clc
clear all
close all



dx=[5150 4000 -1800 -450 -8530 -5270 -450 -1800 4000 5150]/1e0

dy=[0 2000 1180 -185.9 8530*3/100 -5270*3/100 88.1 -1180 -2000 0]/1e0


x=cumsum(dx);
y=cumsum(dy);


x=[0 x];
y=[0 y];


% y=y-1.5;

nodes=[ [1:length(y)].' x.' y.' ];

opt=struct();
opt.ElColor=[0 0 0]
opt.LineWidth=0.75
opt.MarkerSize=NaN;
opt.PlotText=false;
 
% p=inputParser;
% addParameter(p,'MarkerSize',3)
% addParameter(p,'NodeColor',[1 0 0])
% addParameter(p,'ElColor',[0 0 1])
% addParameter(p,'LineWidth',0.75)
% addParameter(p,'PlotText',true,@islogical)
% addParameter(p,'Plotthickness',true,@islogical)
% addParameter(p,'FontSize',6)
% addParameter(p,'shift',false,@islogical)
% addParameter(p,'offset',[],@isnumeric)
% addParameter(p,'mirror',false,@islogical)
% addParameter(p,'hax',[],@ishandle)
% addParameter(p,'interpreter','latex',@ischar)

figure(); 
ha=tight_subplot(1,1,[],[0.1 0.025],[0.1 0.05]);

opt.hax=ha;

scale=100

nodes(:,2:3)=nodes(:,2:3)/1000;

plotthinwalledsection(nodes,[],[],[],[],[],[],opt);

axistight(gca,[0.075 0.3],'x','y');
yticks([0 1 2 3])

% plotscriptmain('h',5,'w',8,'name','FDD','labelsize',6,'ticksize',6,'legendsize',6,'format',{'pdf' 'jpg'});
plotscriptmain('h',3,'w',8,'path','examples/','name','hardanger','labelsize',8,'ticksize',8,'format',{'pdf' 'jpg' 'svg'});

% Thickness=1*ones(length(nodes),1)
% [yc,zc,A_tot,Iy_tot,Iz_tot,Iyz_tot]=SectionParameters(nodes,[],Thickness)
% 
% 
% 
% 
% [A,Cx,Cy,Ix,Iy,Ixy,P]=PolygonMoments([ x.'-yc y.'-zc ])