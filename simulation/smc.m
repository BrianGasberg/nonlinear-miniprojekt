clc
close all
clear all

l = 0.305;
g = -9.81;
Km = 0.0934;
A = 500;
a1 = 8;

dt = 1/200;

a = g/l;

m_p = [0.095 0.105];
alpha_p = [2/2000 5/2000];

c_p = [Km*alpha_p(2)/(m_p(1)*l^2) Km*alpha_p(1)/(m_p(2)*l^2)];
c_hat = mean(c_p);

Kf_p = [0.4 0.6];

b_p = [Kf_p(1)/m_p(1) Kf_p(2)/m_p(1) Kf_p(2)/m_p(2) Kf_p(1)/m_p(2)];

b_hat = mean(b_p);

m = min(m_p);
Kf = max(Kf_p);
alpha = max(alpha_p);

rho1 = [abs(-a/max(c_p) - a/c_hat) abs((a1-max(b_p))/min(c_p) - (a1-b_hat)/c_hat)]
rho2 = (Km*A*(5/2000)/(min(m_p)*l))/min(c_p)

x2_bound = 0.01;
x1_bound = 0.1;
theta = 0.9;
scale = 1;

e = x1_bound/(1-1/(a1*theta))*scale

omega_e = [e/(a1*theta) (1 + 1/(a1*theta))*e];
% omega_e = [x1_bound x2_bound]; 

Rho_wo_dis = rho1(1)*omega_e(1);
beta_0_wo_dis = ceil(Rho_wo_dis) - Rho_wo_dis
beta_wo_dis = beta_0_wo_dis + Rho_wo_dis

Rho = rho1*[sin(omega_e(1)); omega_e(2)] + rho2
beta0 = ceil(Rho) - Rho