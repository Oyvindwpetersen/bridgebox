function [t_a,t_b,I_a,I_b]=LineIntersect(N_a,N_b,n_a,n_b)

%% Calculate intersection of two lines 

% Parametrized vectors:
% v_a=N_a+t_a*n_a
% v_b=N_b+t_b*n_b
% Solve v_a=v_b

% Inputs:
% N_a: coordinates of node 1 
% N_b: coordinates of node 2 
% n_a: direction of vector 1
% n_b: direction of vector 2

% Outputs:
% I_a: coordinates intersection
% I_b: coordinates intersection

% I_a and I_b should be equal

%% Find intersection

A=[-n_a(1) n_b(1) ; -n_a(2) n_b(2) ];
b=[N_a(1)-N_b(1) ; N_a(2)-N_b(2) ];
x=A\b;

t_a=x(1);
t_b=x(2);

I_a=N_a+t_a*n_a;
I_b=N_b+t_b*n_b;

