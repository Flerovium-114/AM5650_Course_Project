function dXdt = bouc_wen_ode(t, X, params, omega, zeta, F0, Omega)

% Inputs:
%   X(1) = x (displacement)
%   X(2) = y (velocity)
%   X(3) = z (hysteretic displacement)

x = X(1);
y = X(2);
z = X(3);

% Hysteresis parameters
beta = params.beta;
gamma = params.gamma;
alpha = params.alpha;
A = params.A;

% Hysteresis term calculation (Eq. 7)
epsilon = sign(y)*sign(z);
dzdt = A*y - (gamma + epsilon*beta)*abs(z)^2*y;

% Main system equations (Eq. 6)
dydt = -2*zeta*omega*y - alpha*omega^2*x - ...
       (1-alpha)*omega^2*z + F0*sin(Omega*t);

dXdt = [y;      % dx/dt
        dydt;   % dy/dt
        dzdt];  % dz/dt
end
