function [yc,zc,A_tot,Iy_tot,Iz_tot,Iyz_tot]=sectionparameters(nodes,elements,thickness)

%% Calculate section properties for twinwalled section
%
% Inputs:
% nodes: matrix with rows [nodenumber, coord_y, coord_z]
% elements: matrix with rows [el_number, nodenumber1, nodenumber2]
% thickness: vector with element thickness corresponding to elements
%
% Outputs:
% yc: center of area
% zc: center of area
% A_tot: area
% Iy_tot: second moment of inertia
% Iz_tot: second moment of inertia
% Iyz_tot: crossed moment of inertia

%% Check inputs

% If elements is empty, assume a closed section with elements from node 1 to node 2
if isempty(elements)
    el_numer=1:size(nodes,1);
    node_number_start=nodes(:,1);
    node_number_end=[ node_number_start(2:end) ; node_number_start(1)];
    elements=[ el_numer.' node_number_start node_number_end];
end

%%

for k=1:size(elements,1)
    
    ind1(k)=find(nodes(:,1)==elements(k,2));
    ind2(k)=find(nodes(:,1)==elements(k,3));
    
    % Vector along  element
    n1{k}=nodes(ind2(k),2:3)-nodes(ind1(k),2:3);
    
    L_el(k)=norm(n1{k});
    
    % Local coordinates
    yc_loc(k)=(nodes(ind1(k),2)+nodes(ind2(k),2))/2;
    zc_loc(k)=(nodes(ind1(k),3)+nodes(ind2(k),3))/2;
    
    %
    A(k)=L_el(k)*thickness(k);
    Sy(k)=zc_loc(k)*A(k);
    Sz(k)=yc_loc(k)*A(k);
    
    
end

zc=sum(Sy)/sum(A);
yc=sum(Sz)/sum(A);

%%

Iy=[];
Iz=[];
Iyz=[];

for k=1:size(elements,1)
    
    % Rotation angle of element relative to y-axis and rotation matrix
    angle_rad(k)=atan2(n1{k}(2),n1{k}(1));
    
    Iy_prime_loc(k)=1/12*thickness(k)^3*L_el(k);
    Iz_prime_loc(k)=1/12*thickness(k)*L_el(k)^3;
    
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

