%% Initialize

init();
clear

p = params_different_gains();
% p = params();
[d, o] = dynamics_and_outputs(p);

%% FL, then MPC on linear system
[X_FL_MPC, T_FL_MPC, U_FL_MPC, X_BAR_FL_MPC, u_FL_MPC, U_FF_FL_MPC] = FL_MPC(p, d, o);

%% MPC-CLF proposed approach
for k = 1:2:15
    p.MPC_CLF.alpha_MPCFL = 0.2*(k-1)/10; % Lf
    p.MPC_CLF.beta_MPCFL = 0.2*(k-1)/10;  % Lg
    [X_Lin_MPC_CLF{k}, T_Lin_MPC_CLF{k}, U_Lin_MPC_CLF{k}, X_K_MPC_CLF{k}, u_Lin_MPC_CLF{k}, U_FF_MPC_CLF{k}] = MPC_Bez(p, d, o, T_FL_MPC, X_FL_MPC, U_FL_MPC);
end


%% Plotting
different_gain_plot();