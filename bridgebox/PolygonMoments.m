function [A,yc,zc,Iy,Iz,Iyz,P]=polygonmoments(coord_y,coord_z)
%% Moments of closed polygon about the origin
%
% Inputs:
% coord_y: y in the counter-clockwise order
% coord_z: z in the counter-clockwise order
%
% Outputs:
% A: area
% yc: center of area 
% zc: center of area
% Iy: second moment of area
% Iz: second moment of area
% Iyz: second crossed moment of area
% P: perimeter circumference
%
%% Moments of area about the origin

coord_yz=[coord_y(:) coord_z(:)]

% At least 3 points required
if max(size(coord_yz))<3
    error('Polygon must have 3 nodes or more');
end

% Ensure coordinates as rows
if size(coord_yz,2)>size(coord_yz,1)
    coord_yz=coord_yz.';
end

% Convert to cell, compute for each cell element
if ~iscell(coord_yz)
    coord_yz_temp=coord_yz;
    clear coord_yz
    coord_yz{1}=coord_yz_temp;
end

for k=1:length(coord_yz)
    [A(k),yc(k),zc(k),Iy(k),Iz(k),Iyz(k),P(k)]=CalculateMoments(coord_yz{k});
end

end

%%
function [A,yc,zc,Iy,Iz,Iyz,P]=CalculateMoments(coord_yz)

% https://en.wikipedia.org/wiki/Second_moment_of_area#Any_polygon

if isempty(coord_yz)
    A=0;
    yc=0;
    zc=0;
    Iy=0;
    Iz=0;
    Iyz=0;
    return
end

y_circ=[coord_yz(:,1) ; coord_yz(1,1)];
z_circ=[coord_yz(:,2) ; coord_yz(1,2)];

temp1=0;
for j=1:length(coord_yz)
	temp1=temp1+y_circ(j)*z_circ(j+1)-y_circ(j+1)*z_circ(j);
end
A=0.5*temp1;
    
temp2=0;
for j=1:length(coord_yz)
	temp2=temp2+(y_circ(j)+y_circ(j+1))*(y_circ(j)*z_circ(j+1)-y_circ(j+1)*z_circ(j));
end
yc=1/(6*A)*temp2;

temp3=0;
for j=1:length(coord_yz)
	temp3=temp3+(z_circ(j)+z_circ(j+1))*(y_circ(j)*z_circ(j+1)-y_circ(j+1)*z_circ(j));
end
zc=1/(6*A)*temp3;

temp4=0;
for j=1:length(coord_yz)
	temp4=temp4+(y_circ(j)*z_circ(j+1)-y_circ(j+1)*z_circ(j)) * (y_circ(j)^2+y_circ(j)*y_circ(j+1)+y_circ(j+1)^2);
end
Iz=1/(12)*temp4;

temp5=0;
for j=1:length(coord_yz)
	temp5=temp5+(y_circ(j)*z_circ(j+1)-y_circ(j+1)*z_circ(j)) * (z_circ(j)^2+z_circ(j)*z_circ(j+1)+z_circ(j+1)^2);
end
Iy=1/(12)*temp5;

temp6=0;
for j=1:length(coord_yz)
	temp6=temp6+(y_circ(j)*z_circ(j+1)-y_circ(j+1)*z_circ(j)) * (y_circ(j)*z_circ(j+1)+2*y_circ(j)*z_circ(j)+2*y_circ(j+1)*z_circ(j+1)+y_circ(j+1)*z_circ(j));
end
Iyz=1/(24)*temp6;

temp7=0;
for j=1:length(coord_yz)
    N1=[y_circ(j) z_circ(j)];
    N2=[y_circ(j+1) z_circ(j+1)];
    v=N2-N1;
	temp7=temp7+norm(v);
end
P=temp7;

end