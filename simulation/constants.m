clc


g = 9.81;
l = 0.305;
kf = 0.5;
m = 0.1;
km = 0.0934;
A = 500;
i = 5/2000*A;

first = km*i

a = g/l;
b = kf/m;

second = (0.4/0.105-b)

c = 1/(l^2*m);

third = (c-(1/(l^2*1.05)))/c


chat = 1/(l^2*0.105);
H = abs((chat-c)*a/c);
epsilon = 0.1;
eta = 1.25;
kappa = 0.05;

A = [ 0 1; 0 -b];
B = [0; c];
[K1, K2] = place(A,B,[-1, -2]);

eig([0 1; -50 -15]);
