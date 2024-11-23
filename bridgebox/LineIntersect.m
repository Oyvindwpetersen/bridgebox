function [t_a,t_b,I_a,I_b]=lineintersect(N_a,N_b,n_a,n_b)
%% Calculate intersection of two lines 
%
% Parametrized vectors:
% v_a=N_a+t_a*n_a
% v_b=N_b+t_b*n_b
% Solve v_a=v_b
%
% Inputs:
% N_a: coordinates of node a
% N_b: coordinates of node b 
% n_a: direction of vector a
% n_b: direction of vector b
%
% Outputs:
% t_a: parameteric value of extension of vector a
% t_b: parameteric value of extension of vector b
% I_a: coordinates intersection
% I_b: coordinates intersection
%
% I_a and I_b should be equal
%
%% Check

if length(N_a)~=2 | length(N_b)~=2 | length(n_a)~=2 | length(n_b)~=2
    N_a
    N_b
    n_a
    n_b
    error('Length of vector must be 2')
end

theta=atan2(norm(cross([n_a 0],[n_b 0])),dot(n_a,n_b));

if theta<0.01
    warning('Angle between vectors is small, almost parallell');
elseif theta<eps
    error('Angle between vectors is zero, no intersection possible');
end    
    
%% Find intersection

A=[-n_a(1) n_b(1) ; -n_a(2) n_b(2) ];
b=[N_a(1)-N_b(1) ; N_a(2)-N_b(2) ];
x=A\b;

t_a=x(1);
t_b=x(2);

I_a=N_a+t_a*n_a;
I_b=N_b+t_b*n_b;