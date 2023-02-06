function [A,Cx,Cy,Ix,Iy,Ixy,P]=PolygonMoments(CoordXY)
%% Moments of closed polygon about the origin

% Inputs:
% CoordXY: matrix with rows [coord_x,coord_y] in the counter-clockwise order
%
% Outputs:
% A: area
% Cx: center of area 
% Cy: center of area
% Ix: second moment of area
% Iy: second moment of area
% Ixy: second crossed moment of area
% P: perimeter circumference
%
%% Moments of area about the origin

% At least 3 points required
if max(size(CoordXY))<3
    error('Polygon must have 3 nodes or more');
end

% Ensure coordinates as rows
if size(CoordXY,2)>size(CoordXY,1)
    CoordXY=CoordXY.';
end

% COnvert to cell, compute for each cell element
if ~iscell(CoordXY)
    CoordXY_temp=CoordXY;
    clear CoordXY
    CoordXY{1}=CoordXY_temp;
end

for k=1:length(CoordXY)
    [A(k),Cx(k),Cy(k),Ix(k),Iy(k),Ixy(k),P(k)]=CalculateMoments(CoordXY{k});
end

end

%%
function [A,Cx,Cy,Ix,Iy,Ixy,P]=CalculateMoments(CoordXY)

% https://en.wikipedia.org/wiki/Second_moment_of_area#Any_polygon

if isempty(CoordXY)
    A=0;
    Cx=0;
    Cy=0;
    Ix=0;
    Iy=0;
    Ixy=0;
    return
end

xCirc=[CoordXY(:,1) ; CoordXY(1,1)];
yCirc=[CoordXY(:,2) ; CoordXY(1,2)];

temp1=0;
for j=1:length(CoordXY)
	temp1=temp1+xCirc(j)*yCirc(j+1)-xCirc(j+1)*yCirc(j);
end
A=0.5*temp1;
    
temp2=0;
for j=1:length(CoordXY)
	temp2=temp2+(xCirc(j)+xCirc(j+1))*(xCirc(j)*yCirc(j+1)-xCirc(j+1)*yCirc(j));
end
Cx=1/(6*A)*temp2;

temp3=0;
for j=1:length(CoordXY)
	temp3=temp3+(yCirc(j)+yCirc(j+1))*(xCirc(j)*yCirc(j+1)-xCirc(j+1)*yCirc(j));
end
Cy=1/(6*A)*temp3;

temp4=0;
for j=1:length(CoordXY)
	temp4=temp4+(xCirc(j)*yCirc(j+1)-xCirc(j+1)*yCirc(j)) * (xCirc(j)^2+xCirc(j)*xCirc(j+1)+xCirc(j+1)^2);
end
Iy=1/(12)*temp4;

temp5=0;
for j=1:length(CoordXY)
	temp5=temp5+(xCirc(j)*yCirc(j+1)-xCirc(j+1)*yCirc(j)) * (yCirc(j)^2+yCirc(j)*yCirc(j+1)+yCirc(j+1)^2);
end
Ix=1/(12)*temp5;

temp6=0;
for j=1:length(CoordXY)
	temp6=temp6+(xCirc(j)*yCirc(j+1)-xCirc(j+1)*yCirc(j)) * (xCirc(j)*yCirc(j+1)+2*xCirc(j)*yCirc(j)+2*xCirc(j+1)*yCirc(j+1)+xCirc(j+1)*yCirc(j));
end
Ixy=1/(24)*temp6;

temp7=0;
for j=1:length(CoordXY)
    N1=[xCirc(j) yCirc(j)];
    N2=[xCirc(j+1) yCirc(j+1)];
    v=N2-N1;
	temp7=temp7+norm(v);
end
P=temp7;

end