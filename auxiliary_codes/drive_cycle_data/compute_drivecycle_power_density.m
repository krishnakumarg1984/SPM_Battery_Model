function [P_density,P_cell,P_batt,P_wheels] = compute_drivecycle_power_density(cycle_name,Mv,Cd,Av,eta_drivetrain,no_of_cells_in_pack,overall_electrochemical_surface_area)
rho = 1.2041; % kg/m^3, At 20°C and 101.325 kPa for dry air (https://en.wikipedia.org/wiki/Density_of_air)
Cr = 0.01;    % Rolling resistance coefficient
g = 9.81;     % m/s^2
Z = 0;        % grade (slope of incline)

load([cycle_name '.mat']); % With sample time of 1 sec
v = drivecycle_data.speed_met_per_sec; % Velocity in m/s
dv_dt = drivecycle_data.acc_met_per_s2;
P_wheels = 0.5*rho*Cd*Av*v.^3 + Cr*Mv*g*v + Mv.*dv_dt.*v + Mv*g*Z*v;
P_batt = P_wheels/eta_drivetrain;
P_cell = P_batt/no_of_cells_in_pack;
P_density = P_cell/overall_electrochemical_surface_area;

end
