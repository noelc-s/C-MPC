function p = params() 

p = struct;

%% Discretization scheme
% p.Disc = 'Exact';
% p.Disc = 'Forward Euler';
% p.Disc = 'Backward Euler';
p.Disc = 'Tustin';

% p.low_level = 'None';
p.low_level = 'CLF';
% p.low_level = 'FL';

%% ODE params
p.ODE.tspan = [0 10];
p.ODE.X0 = [3.14 0];
p.ODE.Xf = [0; 0];

p.sensor_samples_per_sec = 30;
p.noise_mag_sensor = 0.00; % noise added to the dynamics
p.E = 0.001; % size of assumed disturbance invariant

%% Gains
% FL:


p.CLF.eps = 0.9;
p.CLF.alpha_CLF = [p.CLF.eps^2 2*p.CLF.eps];

p.FL.alpha_FL = p.CLF.alpha_CLF;

% CLF Params
p.CLF.F = [0 1; 0 0];
p.CLF.G = [0; 1];
p.CLF.A = p.CLF.F - p.CLF.G*p.CLF.alpha_CLF;
p.CLF.Q = eye(2);
p.CLF.P_lyap = lyap(p.CLF.A', p.CLF.Q);
p.CLF.gamma = min(eig(p.CLF.Q))/max(eig(p.CLF.P_lyap));

% MPC
p.MPC.N = 100;
p.MPC.dt = .1;
p.ll_dt = 0.001;

p.ll_delay = 0; % low_level_state_update delay in ll_ticks
p.xd_delay = floor(0 * p.MPC.dt / p.ll_dt); % desried trajectory update delay in ll_ticks -- this is one mpc timing cycle off

% MPC on Feedback Linearized system:
if p.MPC.N*p.MPC.dt > p.ODE.tspan(2)
    error('There will be no (x_bar, u_bar) because sim was too short!');
end

% MPC-CLF (proposed):
p.MPC_CLF.order = 3;
p.MPC_CLF.alpha_MPCFL = 1; % Lf
p.MPC_CLF.beta_MPCFL = 0;  % Lg
p.MPC_CLF.gamma_MPCFL = 4*max(eig(p.CLF.P_lyap))^3/min(eig(p.CLF.Q))^2; % Robust_level_Set_value
p.MPC_CLF.Gamma_MPCFL = sqrt(p.MPC_CLF.gamma_MPCFL*p.E^2/min(eig(p.CLF.P_lyap))); % E_bar (i.e. sqrt(gamma/lambda_min(P))
p.MPC_CLF.delta_MPCFL = p.MPC_CLF.alpha_MPCFL+sqrt(1+p.MPC_CLF.alpha_MPCFL^2);
p.MPC_CLF.state_stage_cost = 1;
p.MPC_CLF.input_stage_cost = 1;
p.MPC_CLF.aux_stage_cost = 1;

%% Constraints
% State and Input bounds
p.Const.A_in = [1 0; -1 0; 0 1; 0 -1];
p.Const.b_in = [4; 4; 0.5; 0.5];
p.Const.u_max = 20;
p.Const.u_min = -p.Const.u_max;

%% Specifications

% values related to the system dynamics
norm_G = 1; % the norm of g(x) is constant and equal to 1

% values related to the enforced constraints
max_norm_x = norm([max(p.Const.b_in(1),p.Const.b_in(2)) max(p.Const.b_in(3),p.Const.b_in(4))]);

% values related to low level component specification
t_rho_t = @(t) p.Const.u_max*norm_G* t^2 ...
    * exp((p.MPC_CLF.alpha_MPCFL+p.MPC_CLF.beta_MPCFL*p.Const.u_max)*t);
w_mag = @(t) (p.MPC_CLF.alpha_MPCFL + 2*p.MPC_CLF.beta_MPCFL*p.Const.u_max...
    + norm([0 1; -p.FL.alpha_FL]))*t_rho_t(t) + norm_G*p.Const.u_max*t;

%  values related to the timing
T_l_fresh = p.ll_delay * p.ll_dt;
T_m_fresh = p.xd_delay * p.ll_dt;
Delta_T_m = p.MPC.dt - p.ll_dt*(p.MPC.dt-(T_m_fresh + T_l_fresh)) / p.ll_dt;

% mpc related termns
delta_MPC_A_INIT = 0;   % mpc assumes that it gets the true initial state
delta_MPC_G_INIT = p.E; % the desired trajectory will be placed within E of the true inital state

% variation in trajectories
D_x = (p.MPC_CLF.alpha_MPCFL + p.MPC_CLF.beta_MPCFL * p.Const.u_max) * max_norm_x + norm_G * p.Const.u_max;
D_d = p.MPC_CLF.alpha_MPCFL*max_norm_x + p.MPC_CLF.beta_MPCFL * p.Const.u_max;


delta_est_sensor = 0; % the estimator can exactly measure the system state, i.e., \hat x = x

delta_MPC_sensor = delta_est_sensor + T_m_fresh * D_x; % enforced via (4)

delta_MPC_dynamics = 0;

delta_FL_tracking = p.ll_dt*w_mag(p.ll_dt);

delta_FL_dynamics = p.E; % guarantee from MPC

assert(delta_est_sensor +  T_m_fresh * D_x <= delta_MPC_sensor);                  % (4)
assert(delta_MPC_G_INIT + delta_MPC_A_INIT + D_x*T_l_fresh <= delta_FL_dynamics); % (5)
assert(delta_FL_tracking + delta_MPC_dynamics + ...
    (T_m_fresh + T_l_fresh) * D_x + D_d * Delta_T_m <= delta_FL_dynamics);        % (6)
% (7 and 8) are satisfied by definition of terminal set
assert(delta_FL_tracking <= delta_FL_dynamics); % size of robust invariant, must be > w_mag(t_ll);



end