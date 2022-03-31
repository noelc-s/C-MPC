function p = params() 

p = struct;

%% Discretization scheme
p.Disc = 'Exact';
% p.Disc = 'Forward Euler';
% p.Disc = 'Backward Euler';
% p.Disc = 'Tustin';

% p.low_level = 'None';
p.low_level = 'CLF';
% p.low_level = 'FL';

%% ODE params
p.ODE.tspan = [0 15];
p.ODE.X0 = [4 0];
p.ODE.Xf = [0; 0];

p.sensor_samples_per_sec = 10;
p.noise_mag = 0.0;

%% Gains
% FL:
p.FL.alpha_FL = [1 1];

p.CLF.eps = 1.5;
p.CLF.alpha_CLF = [p.CLF.eps^2 2*p.CLF.eps];

% CLF Params
p.CLF.F = [0 1; 0 0];
p.CLF.G = [0; 1];
p.CLF.A = p.CLF.F - p.CLF.G*p.CLF.alpha_CLF;
p.CLF.Q = eye(2);
p.CLF.P_lyap = lyap(p.CLF.A', p.CLF.Q);
p.CLF.gamma = min(eig(p.CLF.Q))/max(eig(p.CLF.P_lyap));

% MPC
p.MPC.N = 10;
p.MPC.dt = 1;
% MPC on Feedback Linearized system:
if p.MPC.N*p.MPC.dt > p.ODE.tspan(2)
    error('There will be no (x_bar, u_bar) because sim was too short!');
end

% MPC-CLF (proposed):
p.MPC_CLF.order = 3;
p.MPC_CLF.alpha_MPCFL = 1; % Lf
p.MPC_CLF.beta_MPCFL = 0;  % Lg
p.MPC_CLF.gamma_MPCFL = 4*max(eig(p.CLF.P_lyap))^3/min(eig(p.CLF.Q))^2; % Robust_level_Set_value
p.MPC_CLF.Gamma_MPCFL = sqrt(p.MPC_CLF.gamma_MPCFL*p.noise_mag^2/min(eig(p.CLF.P_lyap))); % E_bar (i.e. sqrt(gamma/lambda_min(P))
p.MPC_CLF.delta_MPCFL = p.MPC_CLF.alpha_MPCFL+sqrt(1+p.MPC_CLF.alpha_MPCFL^2);
p.MPC_CLF.state_stage_cost = 1;
p.MPC_CLF.input_stage_cost = 1;
p.MPC_CLF.aux_stage_cost = 1;

%% Constraints
% State and Input bounds
p.Const.A_in = [1 0; -1 0; 0 1; 0 -1];
p.Const.b_in = [10; 10; 1; 0.6];
p.Const.u_max = 1.5;
p.Const.u_min = -p.Const.u_max;

end