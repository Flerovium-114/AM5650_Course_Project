function dstatedt = BW_HysteresisFullSystem(t, state, params)
    x = state(1);
    v = state(2);
    z = state(3);

    m = params.m;
    c = params.c;
    k = params.k;
    alpha = params.alpha;
    D = params.D;
    A = params.A;
    beta = params.beta;
    gamma = params.gamma;
    n = params.n;
    F0 = params.F0;
    Omega = params.Omega;
    epsilon = params.epsilon;
    zeta = params.zeta;
    omega = params.omega;

    F = F0 * sin(Omega * t);
    H = alpha * k * x + k * (1 - alpha) * z;
    dxdt = v;
    dvdt = (F - H - c * v) / m;
    %dvdt = -2*zeta*v - alpha*(omega^2)*x - (1 - alpha)*(omega^2)*z;    %%%%
    %dzdt = (1 / D) * (A * v - beta * abs(v) * abs(z)^(n-1) * z - gamma * v * abs(z)^n);
    dzdt = (1/D) * (A - (gamma + epsilon*beta)*(abs(z)^n))*v;   %%%%%
    %dzdt = (1/D) * (A - (0.9)*(abs(z)^n))*v;   % For control

    dstatedt = [dxdt; dvdt; dzdt];
end
