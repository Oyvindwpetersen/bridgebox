function [yc,zc,A_tot,Iy_tot,Iz_tot,Iyz_tot]=SectionParameters(Nodes,Elements,Thickness)

%% Design (closed) bridge section with stiffeners
%
% Inputs:
% NodesBox: [N,2] matrix with rows [Coord_y, Coord_z] for box shape

% Outputs:
% Nodes: [*,3] matrix with rows [NodeNumber Coord_y Coord_z]
% Elements: [*,3] matrix with rows [ElementNumber NodeNumber1 NodeNumber2]
% Thickness: [*,1] matrix with rows [Thickness] corresponding to the elements

%% Parse inputs

% p=inputParser;
% addParameter(p,'Plot',true)


% parse(p,varargin{:});
% DoPlot=p.Results.Plot;


%% Check inputs

% If elements is empty, assume a closed section with elements from node 1 to node 2
if isempty(Elements)
    ElNumber=1:size(Nodes,1);
    NodeNumberStart=Nodes(:,1);
    NodeNumberEnd=[ NodeNumberStart(2:end) ; NodeNumberStart(1)];
    Elements=[ ElNumber.' NodeNumberStart NodeNumberEnd];
end

%%

for k=1:size(Elements,1)
    
    
    ind1(k)=find(Nodes(:,1)==Elements(k,2));
    ind2(k)=find(Nodes(:,1)==Elements(k,3));
    
    % Vector along  element
    n1{k}=Nodes(ind2(k),2:3)-Nodes(ind1(k),2:3);
    
    L_el(k)=norm(n1{k});
    
    % Local coordinates
    yc_loc(k)=(Nodes(ind1(k),2)+Nodes(ind2(k),2))/2;
    zc_loc(k)=(Nodes(ind1(k),3)+Nodes(ind2(k),3))/2;
    
    %
    A(k)=L_el(k)*Thickness(k);
    Sy(k)=zc_loc(k)*A(k);
    Sz(k)=yc_loc(k)*A(k);
    
    
end

zc=sum(Sy)/sum(A);
yc=sum(Sz)/sum(A);

%%

Iy=[];
Iz=[];
Iyz=[];

for k=1:size(Elements,1)
    
    % Rotation angle of element relative to x-axis and rotation matrix
    angle_rad(k)=atan2(n1{k}(2),n1{k}(1));
    
    Iy_prime_loc(k)=1/12*Thickness(k)^3*L_el(k);
    Iz_prime_loc(k)=1/12*Thickness(k)*L_el(k)^3;
    
    Iy_loc(k)=Iy_prime_loc(k)*cos(angle_rad(k))^2+Iz_prime_loc(k)*sin(angle_rad(k))^2;
    Iz_loc(k)=Iy_prime_loc(k)*sin(angle_rad(k))^2+Iz_prime_loc(k)*cos(angle_rad(k))^2;
    Iyz_loc(k)=-(Iy_prime_loc(k)-Iz_prime_loc(k))*sin(angle_rad(k))*cos(angle_rad(k));

    Iy(k)=Iy_loc(k)+(zc_loc(k)-zc)^2*A(k);
    Iz(k)=Iz_loc(k)+(yc_loc(k)-yc)^2*A(k);
    Iyz(k)=Iyz_loc(k)+(yc_loc(k)-yc)*(zc_loc(k)-zc)*A(k);

end

%%


A_tot=sum(A);
Iy_tot=sum(Iy);
Iz_tot=sum(Iz);
Iyz_tot=sum(Iyz);

%
if Iyz_tot<1e-12
    Iyz_tot=0;
end

