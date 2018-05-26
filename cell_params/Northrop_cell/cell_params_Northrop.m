% The struct 'param' is appended with a list of cell-specific parameters

param.I_1C = 60;   % amps (1C-rate current corresponding to cell capacity)

%% Electrode Surface Area in the planar direction
param.A = 2.05268559698939; % m^2; overall active surface area

%% Solid particle radius [m]
param.R_n = 2e-6;
param.R_p = 2e-6;

%% Electrode thicknesses [m]
param.len_n = 88e-6;
param.len_p = 72e-6;

%% Max concentrations of Li in the solid phase [mol/m^3]
param.cs_max_n = 30555;
param.cs_max_p = 51554;

%% Solid diffusion coefficients [m^2 / s]
param.D_n = 3.9e-14;
param.D_p = 1.0e-14;

%% Reaction rate coefficients [m^2.5/(mol^0.5 s)]
param.k_n = 5.031e-11;
param.k_p = 2.334e-11;

%% Electrolyte Conductivity Function
param.ce_init = 1000;   % Electrolyte Li-ions initial concentration [mol/m^3]

%% Stoichiometry Limits
param.theta_max_pos = 0.49550;  % at 100% cell SOC
param.theta_max_neg = 0.85510;  % at 100% cell SOC
param.theta_min_pos = 0.99174;  % at 0% cell SOC
param.theta_min_neg = 0.01429;  % at 0% cell SOC

%% Solid phase volume fraction
param.vol_fraction_neg = 0.4824;
param.vol_fraction_pos = 0.5900;

%% Interfacial surface area Calculations [m^2 / m^3]
param.a_n = 3*param.vol_fraction_neg/param.R_n; % Negative electrode
param.a_p = 3*param.vol_fraction_pos/param.R_p; % Positive electrode

%% OCP function handles
param.compute_Uocp_pos = @compute_Uocp_pos_Northrop;
param.compute_Uocp_neg = @compute_Uocp_neg_Northrop;
