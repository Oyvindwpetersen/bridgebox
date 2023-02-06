function [yc,zc,ys,zs,A,Iyy,Izz,Iyz,It,theta,P]=abdbeam_csprop(E,v,Nodes,Elements,t_mat,Folder,varargin)

%% Function to find cross-sectional properties
% Current functionality: constant thickness, closed section

% Inputs:
% E: E-modulus in N/m^2, can be set to unity (program calculates EA, then divides by E after)
% v: Poisson ratio, can be set to 0.3
% Nodes: [N,2] matrix with y (lateral) and z (vertical) coordinates 
% Elements: [M,3] matrix with rows: [ElementNumber NodeNumber1 NodeNumber2], if empty then closed section  
% t_mat: [M,2] matrix with rows: [ElementNumber thickness]
% Folder: Folder of txt file export, can be set to empty to export to C:\Temp\

% The package https://pypi.org/project/abdbeam/ from python is required

% Outputs:
% yc: coordinate of area center 
% zc: coordinate of area center 
% ys: coordinate of shear center 
% zs: coordinate of shear center 
% A: area
% Iyy: second moment of area about area center 
% Izz: second moment of area about area center 
% Iyz: second moment of area about area center 
% It: shear constant
% theta: angle in deg of principal axes, ccw measured from y-axis
% P: perimeter length

%% Parse inputs

p=inputParser;
addParameter(p,'delete_file',false) % delete txt file from python after analysis
addParameter(p,'plot',true) % plot section 
addParameter(p,'offsetgeometry',false) % offset 
addParameter(p,'offsetdir','in') % offset section outwards or inwards

parse(p,varargin{:});
delete_file=p.Results.delete_file;
do_plot=p.Results.plot;
offsetgeometry=p.Results.offsetgeometry;
offsetdir=p.Results.offsetdir;

%% Process

if size(Nodes,1)<size(Nodes,2)
    Nodes=Nodes.';
end

if length(Nodes)<3; error(''); end

if size(Nodes,2)==2
Nodes=[ [1:length(Nodes)]' Nodes];
end

%%

if isempty(Elements)
    Elements=[ [1:length(Nodes)] ; [1:length(Nodes)] ; [2:length(Nodes) 1]].';
end

if length(t_mat)==1
    t=t_mat;
    t_mat=[ [1:length(Nodes)]' t*ones(length(Nodes),1) ];
end

%%

t_unique=unique(t_mat(:,2));

MaterialId=zeros(length(Elements),1);
% MaterialId(:,1)=Elements(:,1);

for k=1:length(t_unique)
    t_ind{k}=find(t_mat(:,2)==t_unique(k));
    
    ElementNumbers=t_mat(t_ind{k},1);
    for n=1:length(ElementNumbers);
    ind=find(ElementNumbers(n)==Elements(:,1));
    MaterialId(ind,1)=k;
    end
    
end

%% Offset geometry

if offsetgeometry
    
    y=Nodes(:,2); z=Nodes(:,3);    
    y_temp=y-mean(y);
    z_temp=z-mean(z);

    [theta,rho]=cart2pol(y_temp,z_temp);
    theta(theta<0)=theta(theta<0)+2*pi;
    dtheta=diff(theta); 
    if sum(dtheta<0)>sum(dtheta>0)
        direction='cw';
    else
        direction='ccw';
    end

    if strcmpi(offsetdir,'in') & strcmpi(direction,'ccw'); offset_sign=1; end
    if strcmpi(offsetdir,'out') & strcmpi(direction,'ccw'); offset_sign=-1; end

    if strcmpi(offsetdir,'in') & strcmpi(direction,'cw'); offset_sign=-1; end
    if strcmpi(offsetdir,'out') & strcmpi(direction,'cw'); offset_sign=1; end
      

    Nodes_offset=OffsetGeometry(Nodes,offset_sign*t/2);


end


%%

% figure; hold on
% plot(Nodes2(:,1),Nodes2(:,2),'xr')
% plot(Nodes(:,1),Nodes(:,2),'ob')

if isempty(Folder)
    Folder='C:\Temp\';
end

if ~strcmpi(Folder(end),'\')
    Folder=[Folder '\'];
end

% Name of txt file
name_save='abdbeam_csprop_result'

FileNamePython=[Folder 'abdbeam_csprop.py'];


%% Input for python file

data_cell={};

data_cell{end+1}='# Stupid fix for adding path';
data_cell{end+1}='import os';
data_cell{end+1}='base_path=r"C:\Users\oyvinpet\Anaconda3"';
data_cell{end+1}='path=os.pathsep.join([os.path.join(base_path, i) for i in [r"", r"bin", r"Scripts", r"Library\mingw-w64\bin", r"Library\bin"]])';
data_cell{end+1}='os.environ["PATH"]+=os.pathsep+path';

data_cell{end+1}='# ';

data_cell{end+1}='import abdbeam as ab';
data_cell{end+1}='import numpy as np';

data_cell{end+1}=['#Stiffness'];

% data_cell{end+1}=['t=' num2str(t,'%0.3e')];
data_cell{end+1}=['E=float(' num2str(E,'%0.3e') ')'];
data_cell{end+1}=['v=float(' num2str(v,'%0.3e') ')'];

% data_cell{end+1}=['t=' 'float(t)'];
% data_cell{end+1}=['E=' 'float(E)'];
% data_cell{end+1}=['v=' 'float(v)'];
data_cell{end+1}=['G=E/(2*(1+v))'];

data_cell{end+1}='sc=ab.Section()';
data_cell{end+1}='# Create a materials dictionary:';
data_cell{end+1}='mts=dict()';

for k=1:length(t_unique)
    data_cell{end+1}=['t' num2str(k) '=' 'float(' num2str(t_unique(k),'%0.3e') ')'];
    data_cell{end+1}=['mts[' num2str(k) ']=ab.Isotropic('  't' num2str(k) ',E,v)'];
end

data_cell{end+1}='# Create a points dictionary based on Y and Z point coordinates:';
data_cell{end+1}='pts=dict()';

for k=1:length(Nodes)
    y_point=Nodes(k,2);
    z_point=Nodes(k,3);
    data_cell{end+1}=['pts[' num2str(k) ']=ab.Point(' num2str(y_point,'%0.5e') ',' num2str(z_point,'%0.5e') ')'];
end

data_cell{end+1}='# Create a segments dictionary referencing point and material ids:';
data_cell{end+1}='sgs=dict()';

for k=1:length(Elements)

    Node1=Elements(k,2);
    Node2=Elements(k,3);

    NodeIndex1=find(Nodes(:,1)==Node1);
    NodeIndex2=find(Nodes(:,1)==Node2);

    data_cell{end+1}=['# Element ' num2str(Elements(k,1)) ];
    data_cell{end+1}=['sgs[' num2str(k) ']=ab.Segment(' num2str(NodeIndex1) ',' num2str(NodeIndex2) ',' num2str(MaterialId(k,1)) ')'];

    y1=Nodes(NodeIndex1,2);
    z1=Nodes(NodeIndex1,3);

    y2=Nodes(NodeIndex2,2);
    z2=Nodes(NodeIndex2,3);    
    
    ElCoordY{k}=[y1 y2];
    ElCoordZ{k}=[z1 z2];

    L_el=sqrt(diff(ElCoordY{k})^2+diff(ElCoordZ{k})^2);

    if L_el<1e-3
    warning(['Element ' num2str(Elements(k,1)) 'less than 1 mm length' ]);
    end

end

data_cell{end+1}='# Point the dictionaries to the section';
data_cell{end+1}='sc.materials=mts';
data_cell{end+1}='sc.points=pts';
data_cell{end+1}='sc.segments=sgs';

data_cell{end+1}='# Calculate and output section properties';
data_cell{end+1}='sc.calculate_properties()';
data_cell{end+1}='#sc.summary()';

data_cell{end+1}='#ab.plot_section(sc, figsize=(6.4*0.8, 4.8*0.8))';
data_cell{end+1}='#Create a single load case and calculate its internal loads';
data_cell{end+1}='#sc.loads[1]=ab.Load(Vz_s=-100)';
data_cell{end+1}='#sc.calculate_internal_loads()';

data_cell{end+1}='#Plot internal loads';
data_cell{end+1}=['#ab.plot_section_loads(sc, 1, int_load_list=[' '''' 'Nxy' '''' '],title_list=[' '''' 'Abdbeam - Nxy (N/m)' '''' '],figsize=(6.4*0.8, 4.8*0.8))'];


data_cell{end+1}='EA=float(sc.p_c[0,0])';
data_cell{end+1}='EIyy=float(sc.p_c[1,1])';
data_cell{end+1}='EIzz=float(sc.p_c[2,2])';
data_cell{end+1}='EIyz=float(sc.p_c[1,2])';
data_cell{end+1}='GJ=float(sc.p_c[3,3])';

data_cell{end+1}='A=EA/E';
data_cell{end+1}='Iyy=EIyy/E';
data_cell{end+1}='Izz=EIzz/E';
data_cell{end+1}='Iyz=EIyz/E';
data_cell{end+1}='J=GJ/G';

data_cell{end+1}='ValueMatrix=np.array([ sc.yc , sc.zc , sc.ys , sc.zs , A , Iyy, Izz , Iyz, J , sc.principal_axis_angle])';

Folder2=strrep(Folder,'\','\\');
Folder2=strrep(Folder2,'\\\\','\\');
if strcmpi(Folder2(end-1:end),'\\'); Folder2=Folder2(1:end-2); end
    
data_cell{end+1}=['dir_save=' '''' Folder2 ''''];
data_cell{end+1}=['name_save=' '''' name_save ''''];

data_cell{end+1}=['np.savetxt((dir_save+ ''\\'' +name_save+ ''.txt''), ValueMatrix , delimiter='','')'];


%% Write python file

fid_input =fopen(FileNamePython,'wt');
formatSpec ='%s\n';
for k=1:size(data_cell,1)
    fprintf(fid_input,formatSpec,data_cell{k,:});
end
fclose(fid_input);

%% Run python script

% command_str=['python ' Folder 'abdbeam_csprop.py'];

command_str=['cd C:\Users\oyvinpet\Anaconda3\' ' & '  'python ' Folder 'abdbeam_csprop.py'];
[status,commandOut]=system(command_str);

if status~=0
    disp(commandOut);
%     error('Error from python, see above');
end

%% Import

pause(0.01);

% Values=readmatrix([Folder '\' name_save '.txt']);
Values=csvread([Folder '\' name_save '.txt']);

yc=Values(1);
zc=Values(2);

ys=Values(3);
zs=Values(4);

A=Values(5);
Iyy=Values(6);
Izz=Values(7);
Iyz=Values(8);
It=Values(9);
theta=Values(10);

% Round numbers
if abs(yc/sqrt(A))<1e-12
    yc=0;
end

if abs(zc/sqrt(A))<1e-12
    zc=0;
end

if abs(ys/sqrt(A))<1e-12
    ys=0;
end

if abs(zs/sqrt(A))<1e-12
    zs=0;
end

if abs(Iyz/A^2)<1e-12
    Iyz=0;
end

if abs(theta)<1e-3
    theta=0;
end

if delete_file
    delete([Folder '\' name_save '.txt']);
end

if ~do_plot
    return;
end


% for k=1:length(Elements)
%     
%     ind1=Elements(k,2);
%     ind2=Elements(k,3);
%     
%     y1=Nodes(ind1,1);
%     z1=Nodes(ind1,2);
% 
%     y2=Nodes(ind2,1);
%     z2=Nodes(ind2,2);    
%     
%     ElCoordY{k}=[y1 y2];
%     ElCoordZ{k}=[z1 z2];
% 
% end


% Perimeter
for k=1:length(Elements)
    L_el(k)=sqrt(diff(ElCoordY{k})^2+diff(ElCoordZ{k})^2);
end
P=sum(L_el);

%% Plot

figure(); hold on; grid on;
axis image

% y=Nodes(:,1);
% % z=Nodes(:,2);

for k=1:length(Elements)
    
%     ind1=Elements(k,2);
%     ind2=Elements(k,3);
%     
%     y1=Nodes(ind1,2);
%     z1=Nodes(ind1,3);
% 
%     y2=Nodes(ind2,2);
%     z2=Nodes(ind2,3);    
%     
%     ElCoordY{k}=[y1 y2];
%     ElCoordZ{k}=[z1 z2];

    line(ElCoordY{k},ElCoordZ{k},'Color','k','LineWidth',2);
%     text(mean(ElCoordY{k}),mean(ElCoordZ{k}),['EL' num2str(k)],'Color','b','BackGroundColor',[1 1 1],'FontSize',3);
    
end

plot(Nodes(:,2),Nodes(:,3),'ok','MarkerSize',3,'MarkerFaceColor',[0 0 0]);


axistight(gca,[0.25 0.25],'x','y');

% AreaCenterX=sum(a_el_x.*A_el)./sum(A_el);
% AreaCenterZ=sum(a_el_z.*A_el)./sum(A_el);

plot(yc,zc,'xb');
text(yc,zc,'  CA','Color','b','FontSize',8,'BackGroundColor','none');

plot(ys,zs,'xr');
text(ys,zs,'CS  ','Color','r','FontSize',8,'VerticalAlignment','Top','HorizontalAlignment','Right','BackGroundColor','none');

dim_y=max(Nodes(:,2))-min(Nodes(:,2));
dim_z=max(Nodes(:,3))-min(Nodes(:,3));
s=sqrt(dim_y.^2+dim_z.^2)*1.5;

theta_rad=theta*pi/180;
line([ -s/2*cos(theta_rad) s/2*cos(theta_rad)]+yc,[ -s/2*sin(theta_rad) s/2*sin(theta_rad)]+zc,'Color',[0.5 0.5 0.5],'LineStyle','--');

theta2_rad=(theta+90)*pi/180;
line([ -s/2*cos(theta2_rad) s/2*cos(theta2_rad)]+yc,[ -s/2*sin(theta2_rad) s/2*sin(theta2_rad)]+zc,'Color',[0.5 0.5 0.5],'LineStyle','--');

xlabel('y');
ylabel('z');

title_str_cell={};
title_str_cell{end+1}=['A=' num2str(A,'%0.2f')];
title_str_cell{end+1}=['Iyy=' num2str(Iyy,'%0.2f')];
title_str_cell{end+1}=['Izz=' num2str(Izz,'%0.2f')];
title_str_cell{end+1}=['Iyz=' num2str(Iyz,'%0.2f')];
title_str_cell{end+1}=['It=' num2str(It,'%0.2f')];
title_str_cell{end+1}=['theta=' num2str(theta,'%0.2f')];

title_str='';
for k=1:length(title_str_cell)
    title_str=[title_str title_str_cell{k} ', '];
end
title_str=title_str(1:end-2);

title(title_str,'FontWeight','Normal');


end


%% In case of errors, the path to python might not be set
% This can be done every time (or should be fixed permanently)

% open cmd in windows
% run cd C:\Users\oyvinpet\Anaconda3
% run python to check that python exists
% run import os
% run os.environ["PATH"] to check the PATH

% run the code below to add Anaconda to path

% base_path=r"C:\Users\oyvinpet\Anaconda3"
% path=os.pathsep.join([os.path.join(base_path, i) for i in [r"", r"bin", r"Scripts", r"Library\mingw-w64\bin", r"Library\bin"]])
% os.environ["PATH"]+=os.pathsep+path

% https://stackoverflow.com/questions/36778066/importerror-dll-load-failed-when-importing-numpy-installed-in-conda-virtual-env


