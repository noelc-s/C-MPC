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
p.ODE.X0 = [-2 -2];
p.ODE.Xf = [0; 0];

p.sensor_samples_per_sec = 30;
p.noise_mag_sensor = 0.00; % noise added to the dynamics
p.E = 0.05; % size of robust invariant

%% Gains
% FL:
p.FL.alpha_FL = [1 1];

p.CLF.eps = 1;
p.CLF.alpha_CLF = [p.CLF.eps^2 2*p.CLF.eps];

% CLF Params
p.CLF.F = [0 1; 0 0];
p.CLF.G = [0; 1];
p.CLF.A = p.CLF.F - p.CLF.G*p.CLF.alpha_CLF;
p.CLF.Q = eye(2);
p.CLF.P_lyap = lyap(p.CLF.A', p.CLF.Q);
p.CLF.gamma = min(eig(p.CLF.Q))/max(eig(p.CLF.P_lyap));

% MPC
p.MPC.N = 70;
p.MPC.dt = .01;
p.ll_dt = 0.001;

p.ll_delay = 0; % low_level_state_update delay in ll_ticks
p.xd_delay = floor(0 * p.MPC.dt / p.ll_dt); % desried trajectory update delay in ll_ticks -- this is one mpc timing cycle off

% MPC on Feedback Linearized system:
if p.MPC.N*p.MPC.dt > p.ODE.tspan(2)
    error('There will be no (x_bar, u_bar) because sim was too short!');
end

% MPC-CLF (proposed):
p.MPC_CLF.order = 3;
p.MPC_CLF.alpha_MPCFL = .5; % Lf
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
p.Const.b_in = [10; 10; 0.5; 2];
p.Const.u_max = 20;
p.Const.u_min = -p.Const.u_max;

%% Specifications
norm_G = 1;
t_rho_t = @(t) p.Const.u_max*norm_G* t^2 ...
    * exp((p.MPC_CLF.alpha_MPCFL+p.MPC_CLF.beta_MPCFL*p.Const.u_max)*t);
w_mag = @(t) (p.MPC_CLF.alpha_MPCFL + 2*p.MPC_CLF.beta_MPCFL*p.Const.u_max...
    + norm([0 1; -p.FL.alpha_FL]))*t_rho_t(t) + norm_G*p.Const.u_max*t;

T_l_fresh = p.ll_delay * p.ll_dt;
T_m_fresh = p.xd_delay * p.ll_dt;
Delta_T_m = p.MPC.dt - p.ll_dt*(p.ll_dt-(T_m_fresh + T_l_fresh)) / p.ll_dt;

delta_MPC_A_INIT = 0;
delta_MPC_G_INIT = p.E;

max_norm_x = norm([max(p.Const.b_in(1),p.Const.b_in(2)) max(p.Const.b_in(3),p.Const.b_in(4))]);


D_x = (p.MPC_CLF.alpha_MPCFL + p.MPC_CLF.beta_MPCFL * p.Const.u_max) * ...
        max_norm_x + ...
        norm_G * p.Const.u_max;
D_d = p.MPC_CLF.alpha_MPCFL*max_norm_x + p.MPC_CLF.beta_MPCFL * p.Const.u_max;

delta_est_sensor = 0;
delta_MPC_sensor = delta_est_sensor + T_m_fresh * D_x; % from (4)
delta_MPC_dynamics = p.E + p.ll_dt*w_mag(p.ll_dt);

delta_FL_tracking = p.E + p.ll_dt*w_mag(p.ll_dt);

delta_FL_dynamics = p.E; % guarantee from MPC

assert(delta_est_sensor +  T_m_fresh * D_x <= delta_MPC_sensor);                  % (4)
assert(delta_MPC_G_INIT + delta_MPC_A_INIT + D_x*T_l_fresh <= delta_FL_dynamics); % (5)
assert(delta_FL_tracking + delta_MPC_dynamics + ...
    (T_m_fresh + T_l_fresh) * D_x + D_d * Delta_T_m <= delta_FL_dynamics);        % (6)

% (7 and 9) are satisfied by definition of terminal set


assert(p.E >= w_mag(p.ll_dt)); % size of robust invariant, must be > w_mag(t_ll);



end