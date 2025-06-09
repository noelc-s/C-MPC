%% Initialize

init();

p = params();
[d, o] = dynamics_and_outputs(p);

%% FL, then MPC on linear system
[X_FL_MPC, T_FL_MPC, U_FL_MPC, X_BAR_FL_MPC, u_FL_MPC, U_FF_FL_MPC] = FL_MPC(p, d, o);

%% MPC-CLF proposed approach (SOCP)
% p.low_level = 'None';
p.low_level = 'CLF';
[X_Lin_MPC_CLF, T_Lin_MPC_CLF, U_Lin_MPC_CLF, X_K_MPC_CLF, u_Lin_MPC_CLF, U_FF_MPC_CLF] = MPC_Bez(p, d, o, T_FL_MPC, X_FL_MPC, U_FL_MPC);

%% MPC-CLF proposed approach (QP)
% p.low_level = 'None';
p.low_level = 'CLF';
[X_Lin_MPC_CLF_QP, T_Lin_MPC_CLF_QP, U_Lin_MPC_CLF_QP, X_K_MPC_CLF_QP, u_Lin_MPC_CLF_QP, U_FF_MPC_CLF_QP] = MPC_Bez_QP(p, d, o, T_FL_MPC, X_FL_MPC, U_FL_MPC);

