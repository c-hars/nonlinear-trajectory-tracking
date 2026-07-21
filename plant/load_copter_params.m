%% Quadcopter params

qp = struct;
qp.n_rotors = 6;
qp.m = 5.0;
qp.g = 9.81;
qp.l = 0.20;
qp.kF = 0.000015;
qp.kM = qp.kF * 0.01;
qp.I_xx = 0.008;
qp.I_yy = 0.009;
qp.I_zz = 0.015;
qp.Ts = 1/100;

% Motor locations and directions
N = qp.n_rotors;
qp.dirs = (-1).^(1:N);             % -1 = prop rotates CW, 1 = CCW
qp.phi  = (0:N-1)*2*pi/N;
qp.x    =  qp.l * cos(qp.phi + pi/6);  
qp.y    =  qp.l * sin(qp.phi + pi/6);

qp.max_RPM = 10000;
qp.nominal_omegas = compute_omega_bar(qp, ones(1,qp.n_rotors), 'least_squares');
max_omegas = ones(1,qp.n_rotors) * qp.max_RPM*2*pi/60;
max_delta_omegas = max_omegas - qp.nominal_omegas;
qp.max_du = max_delta_omegas;

qp.enabled = [1 1 1 1 1 1]; % all motors working normally