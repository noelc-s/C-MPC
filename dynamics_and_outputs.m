function [d, o] = dynamics_and_outputs(p)
% System Dynamics
syms x1 x2 t
x_sym = [x1 x2];
% f = [x2; x2.^3 + sin(x1)];
f = [x2; sin(x1)];
g = [0; 1];
rng('default')

t_noise = linspace(0,p.ODE.tspan(end),p.sensor_samples_per_sec*p.ODE.tspan(end));
w_noise = p.noise_mag_sensor*2*(rand(1,length(t_noise))-0.5);

% Model of system dynamics to use in controllers
f_model = f;
g_model = g;

% Symbolic gradient for MPC
Df_model = [diff(f,x1) diff(f,x2)];
Dg_model = [diff(g,x1) diff(g,x2)];

% Matlab Function-ify
d.f_func = matlabFunction(f,'Vars',[x_sym]);
d.f_func_w = @(x1,x2,t) d.f_func(x1,x2)+interp1(t_noise, w_noise, t);
d.g_func = matlabFunction(g,'Vars',[x_sym]);
d.g_func_w = @(x1,x2,t) d.g_func(x1,x2)+interp1(t_noise, w_noise, t);
d.f_func_model = matlabFunction(f_model,'Vars',x_sym);
d.g_func_model = matlabFunction(g_model,'Vars',x_sym);
d.Df_func_model = matlabFunction(Df_model,'Vars',x_sym);
d.Dg_func_model = matlabFunction(Dg_model,'Vars',x_sym);

% Define outputs
o.y = x1;
o.Dy = [diff(o.y,x1) diff(o.y,x2)];
o.Lfy = o.Dy*f_model;
o.Lgy = o.Dy*g_model;
o.Lf2y = [diff(o.Lfy,x1) diff(o.Lfy,x2)]*f_model;
o.LgLfy = [diff(o.Lfy,x1) diff(o.Lfy,x2)]*g_model;
o.Lf2y_func = matlabFunction(o.Lf2y,'Vars',x_sym);
o.LgLfy_func = matlabFunction(o.LgLfy,'Vars',x_sym);

% FL
syms v
u = o.LgLfy\(-o.Lf2y + v);
o.FL_u = matlabFunction(u,'Vars',[x1, x2, v]);

% CLF outputs
o.eta = matlabFunction([o.y; o.Lfy],'Vars',x_sym);

end