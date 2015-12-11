
l = 0.305;
g = -9.81;
Km = 0.0934;
A = 500;

dt = 1/200;

a = g/l;

m_p = [0.095 0.105];
alpha_p = [2/2000 5/2000];

c_p = [Km*alpha_p(2)/(m_p(1)*l^2) Km*alpha_p(1)/(m_p(2)*l^2)];
c_hat = mean(c_p);

Kf_p = [0.4 0.6];

b_p = [Kf_p(1)/m_p(1) Kf_p(2)/m_p(1) Kf_p(2)/m_p(2) Kf_p(1)/m_p(2)];

b_hat = mean(b_p);

m = mean(m_p);
Kf = mean(Kf_p);
alpha = mean(alpha_p);


k0_lyap = abs((max(c_p) - c_hat)/c_hat)

theta_lyap = 0.94;

% Assign eigenvalues
k1 = 100;
k2 = 20 - b_hat;

A0 = [0     1;
      -k1 -(k2 + b_hat)];
  
lampda_A0 = eig(A0)

% The Lyapunov functions
P = [3/2 1/2;
    1/2 1/2];

% A0'*P + P*A0

% The pertubation magnitude
% rho1 = 1/c_hat * (abs((a_hat*c - a*c_hat)/c_hat) + abs((c*b_hat - c_hat*b)/c_hat) + k0*sqrt(5));
% rho2 = m_hat;

rho1_lyap = norm((min(c_p) - c_hat)/c_hat^2 *(a - k1)) + norm((b_hat - min(b_p))/c_hat - ((min(c_p) - c_hat)/c_hat^2) * k2)
rho2_lyap = norm(1/c_hat * (Km*A*alpha_p(2))/(min(m_p)*l))

lambda_min = min(eig(P))
lambda_max = max(eig(P))

% The ultimate boundness
b_lyap = 1/2 * sqrt(lambda_max/lambda_min)*sqrt(1/(theta_lyap))

xLimit = norm([0.1 0.01]);

% epsilon
% e_lyap = (xLimit/b_lyap)^2
e_lyap = 0.14