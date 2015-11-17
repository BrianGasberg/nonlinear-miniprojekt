clc
close all
clear all

g = -9.81;
l_min = 0.9;
m_min = 0.5;
k_min = 0.2;


k_hat = 0.1;
l_hat = 1;
m_hat = 0.66;

b_hat = k_hat/m_hat;
a_hat = g/l_hat;
c_hat = 1/(m_hat*l_hat^2);

a = g/l_min;
b = k_min/m_min;
c = 1/(m_min*l_min^2);

k = 0;
l = 0.9;
m = 1.5;

k0 = abs((2.4691 - c_hat)/c_hat);

theta = 0.9;
% k1 = 1;
% k2 = 2 - b_hat;

% Assign eigenvalues
k1 = 40;
k2 = 40 - b_hat;

A0 = [0 1;
      -k1 -(k2 + b_hat)];
  
eig(A0)

% The Lyapunov functions
P = [3/2 1/2;
    1/2 1/2];

% The pertubation magnitude
rho1 = 1/c_hat * (abs((a_hat*c - a*c_hat)/c_hat) + abs((c*b_hat - c_hat*b)/c_hat) + k0*sqrt(5));
rho2 = m_hat;

rho = rho1 + rho2

lambda_min = min(eig(P))
lambda_max = max(eig(P))

% The ultimate boundness
b = sqrt(lambda_max/lambda_min)*sqrt(1/(theta*4))

xLimit = 0.01;

% epsilon
e = (xLimit/b)^2