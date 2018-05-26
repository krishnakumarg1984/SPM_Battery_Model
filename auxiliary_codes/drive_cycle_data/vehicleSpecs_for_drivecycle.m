function [Mv,Cd,Av,eta_drivetrain,no_of_cells_in_pack,overall_electrochemical_surface_area] = vehicleSpecs_for_drivecycle(platform)

pack_mass_chev_bolt = 435.0;   % kg, (not directly used in layer opt calcs) source: Wikipedia and URL: http://media.chevrolet.com/media/us/en/chevrolet/vehicles/bolt-ev/2017.tab1.html
no_of_pax = 2;         % Number of passengers
mass_per_pax = 75.3;   % kg (as per ETA-NTP002 standard, Revision 3, Effective December 1, 2004)
cargo_wt_in_trunk = 0; % kg (weight of cargo, as per above standard)
payload = (no_of_pax*mass_per_pax) + cargo_wt_in_trunk;  % This is the total payload to be carried in the vehicle

Cd = 0.308;            % source: http://www.hybridcars.com/2017-chevy-bolt-ev-is-less-of-a-drag-than-originally-believed/

vehicle_height_inches = 62.9;                    % http://www.chevrolet.com/bolt-ev-electric-vehicle/specs/trims.html
vehicle_height = 0.0254*vehicle_height_inches;   % inches to metres
track_width_inches = 59.1;                       % source:Chevrolet Bolt EV US Website
vehicle_track_width = 0.0254*track_width_inches; % inches to metres
Av = vehicle_height*vehicle_track_width;         % Very conservative estimation, without considering specific body shape/streamlining etc.

eta_drivetrain =  0.75;      % overall drivetrain efficiency (simplified)
cells_in_series = 96;        % no. of series-connected cells in pack

dummy_param{1} = Parameters_init_suppliedSOC_pct(100); % 'dummy' param for lumped mass calculation for Northrop cell
[~,~,mass_Northrop_cell,~] = compute_lumped_mass_and_Cp_avg_for_given_layer_fcn(dummy_param{1}.no_of_layers_Northrop_cell,dummy_param{1});

cells_in_parallel_BEV = 16;       % no. of parallel-connected cells in each module
no_of_cells_in_pack_BEV = cells_in_series*cells_in_parallel_BEV;  % total number of cells in pack
pack_mass_chev_bolt_overhead_BEV = pack_mass_chev_bolt - mass_Northrop_cell*no_of_cells_in_pack_BEV;  % Computing pack mass less the mass of the cells, for the default Bolt configuraiton (kg)

curb_mass_BEV = 1625; % kg ,source: Chevrolet Bolt EV US Website

if strcmp(platform,'BEV')
    curb_mass = curb_mass_BEV;
    Mv = curb_mass + payload + pack_mass_chev_bolt; % Overall vehicle mass (neglecting components with rot. inertia) for all calculations
    no_of_cells_in_pack = no_of_cells_in_pack_BEV;  % For function return
    
elseif strcmp(platform,'PHEV')
    engine_mass = 98; % 98 kg Ecotec engine mass added for PHEV
    curb_mass = (curb_mass_BEV + engine_mass); % kg,
    cells_in_parallel_PHEV = 1;  % no. of parallel-connected cells in pack
    no_of_cells_in_pack_PHEV = cells_in_series*cells_in_parallel_PHEV;  % total number of cells in pack
    no_of_cells_in_pack = no_of_cells_in_pack_PHEV; % For function return
    pack_mass_chev_bolt_overhead_PHEV = pack_mass_chev_bolt_overhead_BEV*(cells_in_parallel_PHEV/cells_in_parallel_BEV);  % Computing pack mass less the mass of the cells (kg)
    Mv = curb_mass + payload + pack_mass_chev_bolt_overhead_PHEV + (mass_Northrop_cell*no_of_cells_in_pack); % Overall vehicle mass (neglecting components with rot. inertia) for all calculations
else
    error('\nVehicle type must be either BEV or PHEV !....\n');
end

overall_electrochemical_surface_area = dummy_param{1}.overall_surface_area_for_given_layers;
end