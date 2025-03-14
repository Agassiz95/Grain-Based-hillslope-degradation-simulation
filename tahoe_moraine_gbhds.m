%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Tahoe Moraine Model                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is the best we will get without fudging the numbers!
tStart = tic;

% Duration of simulation
maxtime = 45000; % Years (Tahoe: 45000; Mono Basin: 100000)

% Setting up initial hillslope profile
array_length = 6000;
xdistance = zeros(1, array_length);
ydistance = zeros(1, array_length);
%alfa = 28; % Average angle for young moraines via Putkonen & O'Neal, 2006
alfa = 30;

% Flat part at the base (left side)
for counter = 1:1900 % originally 1800
    xdistance(counter) = counter;
    ydistance(counter) = 0;  % Flat part at the base
end

% Increasing slope for triangular shape
for counter = 1901:3000
    xdistance(counter) = counter;
    ydistance(counter) = (counter - 1900) * tan((alfa / 360) * 2 * pi);  % Linear increase in height
end

% Decreasing slope (mirror image)
for counter = 3001:4100
    xdistance(counter) = counter;
    ydistance(counter) = ydistance(6001 - counter);  % Symmetry in the triangle
end

% Flat part at the base (right side)
for counter = 4101:array_length % 4201
    xdistance(counter) = counter;
    ydistance(counter) = 0;  % Back to flat base
end

slope_original = ydistance;
node_count = length(xdistance);

%%%%% Difficult to find values doe to unexpected measured grain size
%%%%% distribution. Latin Hypercube Sampling was used to getan idea of
%%%%% what the best starting solution would may look like.

% Transport factor values --> Real values are halved since this only
% simulates one side of the hillslope!
KSmall = 0.0125;   % Best: 0.001 --> RMSE = 2.3
KMedium = 0.005;  % Best: 0.001 --> RMSE = 6.1
KMedLar = 0.005;  % Best: 0.001 --> RMSE = 5.0
KLarge = 0.0025;  % Best: 0.0009 --> RMSE = 3.5
KVeryLarge = 0.002; % Best: 0.0006 --> RMSE = 8.2
KTotal = KSmall + KMedium + KMedLar + KLarge + KVeryLarge; % Best: RMSE = 5.0, KTotal value = 0.0045

% Grain Transport model; Best RMSE = 1.8604 --> uses above best K-values and below best Wx values
K_nonlocal = 0.014; % Best: 0.0026 --> RMSE = 1.8734
K_diffuse = 0.014; % Best: 0.0026 --> RMSE = 1.8285


% Comminution rates for each grain size (%/half year) Tim had 0.00001, so each year 0.002% of the particles become smaller
% Putkonen et al. 2013 says boulders break apart 250x faster than the
% broken components. This is in antarctic dry valleys though So the rate may be
% higher than this. 
% To find the break down rate adjust the largest grain size break down rate
% such that by the end of the model run the grain size distribution near
% the crest matches the observed values.

very_Large_wr = 0.00004; % Best: 0.00004 || 
large_wr = 0.000002; % Best: 0.000002 ||
medlar_wr = 0.000001; % Best: 0.000002 || 
medium_wr = 0.000001; % Best: 0.000002 || 
small_wr = 0.0000; % Best: 0.0000

% The comminution rates are way too high on the left hand values. The K-values associated with this
% the comminution rates hardly lower the crest of the moraine! ie. KTotal
% is way too low!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Grain transport & distribution diffusion model           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time increment
dt = 0.1; % Years normally set to 0.5

% Horizontal increment
dx = 0.1; % meters

% Active volume
active_volume = 1; % decimeters

% Q-values
qSmall = zeros(1, array_length);
qMedium = zeros(1, array_length);
qMedLar = zeros(1, array_length);
qLarge = zeros(1, array_length);
qVeryLarge = zeros(1, array_length);

% Initial grain percents- Taken from tim's data for the Tahoe Moraine 
% Crest, and adjusted by 5% (from the fines) to simulate conditions prior
% to erosion and comminution the assumption is that there is little/no 

initial_percent_small = 0.4; 
initial_percent_medium = 0.13;
initial_percent_medlar = 0.13;
initial_percent_large = 0.14; 
initial_percent_Very_Large = 0.2;

% Initialize active percent of each grain size
active_percent_small = ones(1, array_length) * initial_percent_small;
active_percent_medium = ones(1, array_length) * initial_percent_medium;
active_percent_medlar = ones(1, array_length) * initial_percent_medlar;
active_percent_large = ones(1, array_length) * initial_percent_large;
active_percent_Very_Large = ones(1, array_length) * initial_percent_Very_Large;

% Calculate the number of time steps when rem(time,10) == 0
num_appends = floor(maxtime / 10) + 1; % Adding 1 to include the initial time (time = 0)

% Pre-allocate arrays to store weight percent over time
small_footslope_array = zeros(1, num_appends);
medium_footslope_array = zeros(1, num_appends);
medlar_footslope_array = zeros(1, num_appends);
large_footslope_array = zeros(1, num_appends);
very_Large_footslope_array = zeros(1, num_appends);

small_midslope_array = zeros(1, num_appends);
medium_midslope_array = zeros(1, num_appends);
medlar_midslope_array = zeros(1, num_appends);
large_midslope_array = zeros(1, num_appends);
very_Large_midslope_array = zeros(1, num_appends);

small_crest_array = zeros(1, num_appends);
medium_crest_array = zeros(1, num_appends);
medlar_crest_array = zeros(1, num_appends);
large_crest_array = zeros(1, num_appends);
very_Large_crest_array = zeros(1, num_appends);

time_array = zeros(1, num_appends);
time_array_rates = zeros(1, num_appends);

save_interval = 10; % Save the profile every 10 years
save_idx = 1; % Index for storing profiles

% Initialize index for appending data
idx = 1;
idx_rates = 1;

% Imports raw profile to now use
Grain_transport_profile = slope_original;

dydx = zeros(1, array_length);
for i = 2:(array_length - 1)
    dydx(i) = (Grain_transport_profile(i + 1) - Grain_transport_profile(i - 1)) / (2 * dx);
end

fprintf('Begin Grain Transport and distribution Simulation\n');

% Initialize arrays to store time and crest height
num_steps = floor(maxtime / save_interval) + 1; % Number of intervals to save the profile
crest_height_over_time = zeros(1, num_steps);    % Stores crest height at each interval
time_points = 0:save_interval:maxtime;           % Corresponding time points



%%%%%%%%%%%%%%%%%%%%

% Main simulation loop
for time = 0:dt:maxtime

    % Display time every 1000 years
    if rem(time, 1000) == 0
        fprintf('Year: %d\n', time);
    end

    % Calculate hillslope gradients for all segments
    for i = 2:array_length
        if Grain_transport_profile(i) - Grain_transport_profile(i - 1) > 1e-7
            dydx(i) = (Grain_transport_profile(i) - Grain_transport_profile(i - 1)) / dx;
        else
            dydx(i) = 0;
        end
    end

    % Calculate slope angle in degrees
    slope_angle = abs(-atand(dydx / 10));

    % Modify initial percentages due to weathering
    initial_percent_large = initial_percent_large * (1 - large_wr);
    initial_percent_small = initial_percent_small * (1 + large_wr);

    % Call the function for each grain size
    qSmall = calculate_flux(KSmall, dydx, dt, dx);
    qMedium = calculate_flux(KMedium, dydx, dt, dx);
    qMedLar = calculate_flux(KMedLar, dydx, dt, dx);
    qLarge = calculate_flux(KLarge, dydx, dt, dx);
    qVeryLarge = calculate_flux(KVeryLarge, dydx, dt, dx);
    qTotal = calculate_flux(KTotal, dydx, dt, dx);

    % Initialize added and lost volume arrays
    added_volume_small = zeros(1, array_length);
    added_volume_medium = zeros(1, array_length);
    added_volume_medlar = zeros(1, array_length);
    added_volume_large = zeros(1, array_length);
    added_volume_Very_Large = zeros(1, array_length);

    lost_volume_small = zeros(1, array_length);
    lost_volume_medium = zeros(1, array_length);
    lost_volume_medlar = zeros(1, array_length);
    lost_volume_large = zeros(1, array_length);
    lost_volume_Very_Large = zeros(1, array_length);

    % Factor in surface availability
    qzSmall = zeros(1, array_length);
    qzMedium = zeros(1, array_length);
    qzMedLar = zeros(1, array_length);
    qzLarge = zeros(1, array_length);
    qzVery_Large = zeros(1, array_length);
    qzTotal = zeros(1, array_length);

   % Call the flux calculation function
    [qxSmall, qxMedium, qxMedLar, qxLarge, qxVeryLarge, qxTotal] = calculate_fluxes(array_length, ...
        qTotal, qSmall, qMedium, qMedLar, qLarge, qVeryLarge, ...
        active_percent_small, active_percent_medium, active_percent_medlar, active_percent_large, active_percent_Very_Large);

    % Grain transport for a given hillslope angle
    [added_volume_small, lost_volume_small] = calculate_grain_movement(slope_angle, qxSmall, added_volume_small, lost_volume_small);
    [added_volume_medium, lost_volume_medium] = calculate_grain_movement(slope_angle, qxMedium, added_volume_medium, lost_volume_medium);
    [added_volume_medlar, lost_volume_medlar] = calculate_grain_movement(slope_angle, qxMedLar, added_volume_medlar, lost_volume_medlar);
    [added_volume_large, lost_volume_large] = calculate_grain_movement(slope_angle, qxLarge, added_volume_large, lost_volume_large);
    [added_volume_Very_Large, lost_volume_Very_Large] = calculate_grain_movement(slope_angle, qxVeryLarge, added_volume_Very_Large, lost_volume_Very_Large);


    % Net volume calculations
    added_volume_total = added_volume_small + added_volume_medium + added_volume_medlar + added_volume_large + added_volume_Very_Large;
    lost_volume_total = lost_volume_small + lost_volume_medium + lost_volume_medlar + lost_volume_large + lost_volume_Very_Large;
    net_volume_small = added_volume_small + lost_volume_small;
    net_volume_medium = added_volume_medium + lost_volume_medium;
    net_volume_medlar = added_volume_medlar + lost_volume_medlar;
    net_volume_large = added_volume_large + lost_volume_large;
    net_volume_Very_Large = added_volume_Very_Large + lost_volume_Very_Large;
    net_volume_total = net_volume_small + net_volume_medium + net_volume_medlar + net_volume_large + net_volume_Very_Large;

    % Update the slope profile
    Grain_transport_profile(1:array_length) = Grain_transport_profile(1:array_length) + net_volume_total(1:array_length);

    % Accumulation and erosion effects
    new_volume_small = zeros(1, array_length);
    new_volume_medium = zeros(1, array_length);
    new_volume_medlar = zeros(1, array_length);
    new_volume_large = zeros(1, array_length);
    new_volume_Very_Large = zeros(1, array_length);
    new_volume_total = zeros(1, array_length);

    new_w_volume_Very_Large = zeros(1, array_length);
    new_w_volume_large = zeros(1, array_length);
    new_w_volume_medlar = zeros(1, array_length);
    new_w_volume_medium = zeros(1, array_length);
    new_w_volume_small = zeros(1, array_length);
    new_w_volume_total = zeros(1, array_length);

    for i = 5:array_length
        % Accumulation
        if net_volume_total(i) < active_volume && net_volume_total(i) > 0
            new_volume_small(i) = ((active_volume - net_volume_total(i)) * active_percent_small(i)) + net_volume_small(i);
            new_volume_medium(i) = ((active_volume - net_volume_total(i)) * active_percent_medium(i)) + net_volume_medium(i);
            new_volume_medlar(i) = ((active_volume - net_volume_total(i)) * active_percent_medlar(i)) + net_volume_medlar(i);
            new_volume_large(i) = ((active_volume - net_volume_total(i)) * active_percent_large(i)) + net_volume_large(i);
            new_volume_Very_Large(i) = ((active_volume - net_volume_total(i)) * active_percent_Very_Large(i)) + net_volume_Very_Large(i);

            % Normalize the volumes to ensure the fractions add up to 1
            new_volume_total(i) = new_volume_small(i) + new_volume_medium(i) + new_volume_medlar(i) + new_volume_large(i) + new_volume_Very_Large(i);
            if new_volume_total(i) > 0
                active_percent_small(i) = new_volume_small(i) / new_volume_total(i);
                active_percent_medium(i) = new_volume_medium(i) / new_volume_total(i);
                active_percent_medlar(i) = new_volume_medlar(i) / new_volume_total(i);
                active_percent_large(i) = new_volume_large(i) / new_volume_total(i);
                active_percent_Very_Large(i) = new_volume_Very_Large(i) / new_volume_total(i);
            end

            % Add assertion
            assert(abs(active_percent_small(i) + active_percent_medium(i) + ...
                active_percent_medlar(i) + active_percent_large(i) + active_percent_Very_Large(i) - 1) < 1e-6, ...
                'Error: The sum of the percentages in accumulation does not equal 1 at index %d', i);

            % Weathering
            new_w_volume_Very_Large(i) = new_volume_Very_Large(i) * (1 - very_Large_wr);
            new_w_volume_large(i) = new_volume_large(i) + (very_Large_wr * new_volume_Very_Large(i));
            new_w_volume_medlar(i) = new_volume_medlar(i) + (large_wr * new_volume_large(i));
            new_w_volume_medium(i) = new_volume_medium(i) + (medlar_wr * new_w_volume_medlar(i));
            new_w_volume_small(i) = new_volume_small(i) + (medium_wr * new_w_volume_medium(i));
            new_w_volume_total(i) = new_w_volume_small(i) + new_w_volume_medium(i) + new_w_volume_medlar(i) + new_w_volume_large(i) + new_w_volume_Very_Large(i);
            
            % Normalizing weathering effects
            if new_w_volume_total(i) > 0
                active_percent_small(i) = new_w_volume_small(i) / new_w_volume_total(i);
                active_percent_medium(i) = new_w_volume_medium(i) / new_w_volume_total(i);
                active_percent_medlar(i) = new_w_volume_medlar(i) / new_w_volume_total(i);
                active_percent_large(i) = new_w_volume_large(i) / new_w_volume_total(i);
                active_percent_Very_Large(i) = new_w_volume_Very_Large(i) / new_w_volume_total(i);
            end

            % Add assertion
            assert(abs(active_percent_small(i) + active_percent_medium(i) + ...
                active_percent_medlar(i) + active_percent_large(i) + active_percent_Very_Large(i) - 1) < 1e-6, ...
                'Error: The sum of the percentages in weathering does not equal 1 at index %d', i);

        end

        % Erosion
        if net_volume_total(i) < 0
            new_volume_small(i) = (active_volume * active_percent_small(i)) + net_volume_small(i) - (net_volume_total(i) * initial_percent_small);
            new_volume_medium(i) = (active_volume * active_percent_medium(i)) + net_volume_medium(i) - (net_volume_total(i) * initial_percent_medium);
            new_volume_medlar(i) = (active_volume * active_percent_medlar(i)) + net_volume_medlar(i) - (net_volume_total(i) * initial_percent_medlar);
            new_volume_large(i) = (active_volume * active_percent_large(i)) + net_volume_large(i) - (net_volume_total(i) * initial_percent_large);
            new_volume_Very_Large(i) = (active_volume * active_percent_Very_Large(i)) + net_volume_Very_Large(i) - (net_volume_total(i) * initial_percent_Very_Large);

            % Normalize the volumes to ensure the fractions add up to 1
            new_volume_total(i) = new_volume_small(i) + new_volume_medium(i) + new_volume_medlar(i) + new_volume_large(i) + new_volume_Very_Large(i);
            if new_volume_total(i) > 0
                active_percent_small(i) = new_volume_small(i) / new_volume_total(i);
                active_percent_medium(i) = new_volume_medium(i) / new_volume_total(i);
                active_percent_medlar(i) = new_volume_medlar(i) / new_volume_total(i);
                active_percent_large(i) = new_volume_large(i) / new_volume_total(i);
                active_percent_Very_Large(i) = new_volume_Very_Large(i) / new_volume_total(i);
            end


            % Add assertion
            assert(abs(active_percent_small(i) + active_percent_medium(i) + ...
                active_percent_medlar(i) + active_percent_large(i) + active_percent_Very_Large(i) - 1) < 1e-6, ...
                'Error: The sum of the percentages in erosion does not equal 1 at index %d', i);

            
            % Weathering
            new_w_volume_Very_Large(i) = new_volume_Very_Large(i) * (1 - very_Large_wr);
            new_w_volume_large(i) = new_volume_large(i) + (very_Large_wr * new_volume_Very_Large(i));
            new_w_volume_medlar(i) = new_volume_medlar(i) + (large_wr * new_volume_large(i));
            new_w_volume_medium(i) = new_volume_medium(i) + (medlar_wr * new_w_volume_medlar(i));
            new_w_volume_small(i) = new_volume_small(i) + (medium_wr * new_w_volume_medium(i));
            new_w_volume_total(i) = new_w_volume_small(i) + new_w_volume_medium(i) + new_w_volume_medlar(i) + new_w_volume_large(i) + new_w_volume_Very_Large(i);


            active_percent_small(i) = new_w_volume_small(i) / new_w_volume_total(i);
            active_percent_medium(i) = new_w_volume_medium(i) / new_w_volume_total(i);
            active_percent_medlar(i) = new_w_volume_medlar(i) / new_w_volume_total(i);
            active_percent_large(i) = new_w_volume_large(i) / new_w_volume_total(i);
            active_percent_Very_Large(i) = new_w_volume_Very_Large(i) / new_w_volume_total(i);
        end
    end

    % Store data at specified intervals
    if rem(time, 10) == 0
        % Calculate weight percentages at specific positions
        small_footslope = active_percent_small(1670) * 100; % Footslope = 1860 for Tahoe, 1570
        medium_footslope = active_percent_medium(1670) * 100;
        medlar_footslope = active_percent_medlar(1670) * 100;
        large_footslope = active_percent_large(1670) * 100;
        very_Large_footslope = active_percent_Very_Large(1670) * 100;

        small_midslope = active_percent_small(2400) * 100; % Midslope = 2400 for Tahoe
        medium_midslope = active_percent_medium(2400) * 100;
        medlar_midslope = active_percent_medlar(2400) * 100;
        large_midslope = active_percent_large(2400) * 100;
        very_Large_midslope = active_percent_Very_Large(2400) * 100;

        small_crest = active_percent_small(2940) * 100; % crest = 2950 for Tahoe
        medium_crest = active_percent_medium(2940) * 100;
        medlar_crest = active_percent_medlar(2940) * 100;
        large_crest = active_percent_large(2940) * 100;
        very_Large_crest = active_percent_Very_Large(2940) * 100;

        % Store values in pre-allocated arrays
        small_footslope_array(idx) = small_footslope;
        medium_footslope_array(idx) = medium_footslope;
        medlar_footslope_array(idx) = medlar_footslope;
        large_footslope_array(idx) = large_footslope;
        very_Large_footslope_array(idx) = very_Large_footslope;

        small_midslope_array(idx) = small_midslope;
        medium_midslope_array(idx) = medium_midslope;
        medlar_midslope_array(idx) = medlar_midslope;
        large_midslope_array(idx) = large_midslope;
        very_Large_midslope_array(idx) = very_Large_midslope;

        small_crest_array(idx) = small_crest;
        medium_crest_array(idx) = medium_crest;
        medlar_crest_array(idx) = medlar_crest;
        large_crest_array(idx) = large_crest;
        very_Large_crest_array(idx) = very_Large_crest;

        time_array(idx) = time;

        % Increment index
        idx = idx + 1;
    end

    mid_idx = floor(array_length / 2);
    for i = 1:mid_idx
        symmetric_value = (Grain_transport_profile(i) + Grain_transport_profile(array_length - i + 1)) / 2;
        Grain_transport_profile(i) = symmetric_value;
        Grain_transport_profile(array_length - i + 1) = symmetric_value;
    end

   % Store the crest height at intervals
    if rem(time, save_interval) == 0
        crest_height_over_time(save_idx) = Grain_transport_profile(3000); % Save height at crest index
        save_idx = save_idx + 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Indices for midslope, crest, and footslope
    footslope_idx = 1570;
    midslope_idx = 2400;
    crest_idx = 2940;
    
    % Store transport rates over time
    if rem(time, 10) == 0  % Store every 10 years

        time_array_rates(idx_rates) = time;

        small_footslope_rate(idx_rates) = qSmall(footslope_idx);
        small_midslope_rate(idx_rates) = qSmall(midslope_idx);
        small_crest_rate(idx_rates) = qSmall(crest_idx);
    
        medium_footslope_rate(idx_rates) = qMedium(footslope_idx);
        medium_midslope_rate(idx_rates) = qMedium(midslope_idx);
        medium_crest_rate(idx_rates) = qMedium(crest_idx);
    
        medlar_footslope_rate(idx_rates) = qMedLar(footslope_idx);
        medlar_midslope_rate(idx_rates) = qMedLar(midslope_idx);
        medlar_crest_rate(idx_rates) = qMedLar(crest_idx);
    
        large_footslope_rate(idx_rates) = qLarge(footslope_idx);
        large_midslope_rate(idx_rates) = qLarge(midslope_idx);
        large_crest_rate(idx_rates) = qLarge(crest_idx);
    
        veryLarge_footslope_rate(idx_rates) = qVeryLarge(footslope_idx);
        veryLarge_midslope_rate(idx_rates) = qVeryLarge(midslope_idx);
        veryLarge_crest_rate(idx_rates) = qVeryLarge(crest_idx);
        
        idx_rates = idx_rates+1;
    end





end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic diffusion profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time increment
dt = 0.2; % Years, normally set to 0.5

% Horizontal increment
dx = 0.1; % meters

slope_height = slope_original;  % Initial height profile
time_steps = maxtime / dt + 1;  % Total number of time steps

% Matrix to store slope height at each time step
slope_evolution = zeros(time_steps, length(slope_height));

% Save the initial slope
slope_evolution(1, :) = slope_height;

fprintf('Begin The Linear Diffusion Simulation. \n')

% Diffusion process
for t = 2:time_steps
    % Compute the second-order difference of slope_height
    slope_diff = diff(slope_height, 2);  % Central difference

    % Update slope_height for nodes 2 to node_count-1
    slope_height(2:end - 1) = slope_height(2:end - 1) - (-K_diffuse * (slope_diff * dt / dx^2));


     % Store the updated slope profile
    slope_evolution(t, :) = slope_height;

end

fprintf('The Linear Diffusion Simulation Has Ended. \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Crest position for linear diffusion model
[~, crest_position_diffusion] = max(slope_evolution(1, :));
crest_elevation_diffusion = slope_evolution(:, crest_position_diffusion);
time_diffusion = (0:time_steps-1) * dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Start Non-Local model                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for the non-local model
dx = 1;     % Spatial step (m)
dt = 10;    % Time step (years), usually set to 10
steps = maxtime / dt;  % Number of time steps


% Initial slope profile
z0_nonlocal = slope_original;  % Same initial slope as the diffusion model
z_nonlocal = z0_nonlocal;  % Non-local model elevation profile

% Weighting function for non-local influence (power law decay)
alpha = 2;  % Exponent for power-law decay (adjust as needed)
L = 50;  % Characteristic length scale

omega = @(x_prime) 1 ./ (1 + abs(x_prime) / L) .^ alpha;

% Precompute distances and weights for non-local model
x_matrix_nonlocal = repmat(xdistance', 1, length(xdistance));

distances_nonlocal = abs(x_matrix_nonlocal - x_matrix_nonlocal');  % Calculate distances between points
weights_nonlocal = sparse(omega(abs(distances_nonlocal)));


fprintf('Begin Non-Local Simulation:\n')

% Matrix to store slope profiles at each time step
z_evolution_nonlocal = zeros(steps, length(z_nonlocal));

% Save the initial slope
z_evolution_nonlocal(1, :) = z_nonlocal;

% Time-stepping loop for non-local model
for t = 2:steps
    % Non-local model
    q_nonlocal = sediment_flux(z_nonlocal, weights_nonlocal, dx, K_nonlocal);  % Use precomputed weights and KTotal

    dzdt_nonlocal = -gradient(q_nonlocal, dx);  % Rate of change of elevation

    z_nonlocal = z_nonlocal + dzdt_nonlocal * dt;  % Update elevation

    % Store the updated slope profile
    z_evolution_nonlocal(t, :) = z_nonlocal;

    % Print progress every 1000 years (every 100 time steps)
    if mod(t, 100) == 0
        fprintf('Time step: %d, Time: %d years\n', t, t * dt);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the crest position (peak of the initial profile)
[~, crest_position] = max(z_evolution_nonlocal(1, :));

% Extract elevation values at the crest position over time
crest_elevation_nonlocal = z_evolution_nonlocal(:, crest_position);

% Define the time vector for plotting
time_nonlocal = (0:steps-1) * dt;  % in years

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         End Simulation                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% End timing of simulation 
tEnd = toc(tStart)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          Begin Grain Size Distribution Error Calculations           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert active percent to out of 100 for aesthetics
active_percent_small = active_percent_small * 100;
active_percent_medium = active_percent_medium * 100;
active_percent_medlar = active_percent_medlar * 100;
active_percent_large = active_percent_large * 100;
active_percent_Very_Large = active_percent_Very_Large * 100;

% Main script for RMSE calculations
% Tahoe moraine data points (Foot, Mid, Crest) for each grain size
points = [1670, 2400, 2940]; % Location based on % distance from crest whereby the bulk sample grain size distributions were calculated
% 1570

% Observed values for Tahoe moraine 
observed_small = [active_percent_small(1670), active_percent_small(2400), active_percent_small(2940)];
observed_medium = [active_percent_medium(1670), active_percent_medium(2400), active_percent_medium(2940)];
observed_medlar = [active_percent_medlar(1670), active_percent_medlar(2400), active_percent_medlar(2940)];
observed_large = [active_percent_large(1670), active_percent_large(2400), active_percent_large(2940)];
observed_Very_Large = [active_percent_Very_Large(1670), active_percent_Very_Large(2400), active_percent_Very_Large(2940)];

% Expected values for Tahoe moraine (XSmall, Small, Medium, MedLar, Large)
% [footslope, midslope, crest]
expected_small = [42, 37, 38];
expected_medium = [30, 24, 5];
expected_medlar = [17, 15, 3];
expected_large = [10, 8, 1];
expected_Very_Large = [1, 16, 53];

% Calculate RMSE for each grain size
RMSE_Small_Young = calculate_rmse(expected_small, observed_small);
RMSE_Medium_Young = calculate_rmse(expected_medium, observed_medium);
RMSE_MedLar_Young = calculate_rmse(expected_medlar, observed_medlar);
RMSE_Large_Young = calculate_rmse(expected_large, observed_large);
RMSE_Very_Large_Young = calculate_rmse(expected_Very_Large, observed_Very_Large);

% Total RMSE for Tahoe moraine
RMSE_Total_Young = (RMSE_Small_Young + RMSE_Medium_Young + ...
                        RMSE_MedLar_Young + RMSE_Large_Young + RMSE_Very_Large_Young) / 5;

% Display the RMSE results using fprintf for formatted output
fprintf('Root Mean Square Error (RMSE) for GSD on Tahoe moraine:\n');
fprintf('RMSE_Small_Young: %.1f\n', RMSE_Small_Young);
fprintf('RMSE_Medium_Young: %.1f\n', RMSE_Medium_Young);
fprintf('RMSE_MedLar_Young: %.1f\n', RMSE_MedLar_Young);
fprintf('RMSE_Large_Young: %.1f\n', RMSE_Large_Young);
fprintf('RMSE_Very_Large_Young: %.1f\n', RMSE_Very_Large_Young);
fprintf('RMSE_Total_Young: %.1f\n', RMSE_Total_Young);

fprintf('Ending Grain Size Distribution Error Calculation\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Begin Hillslope Profile Error Calculations               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load the interpolated X values from the Tahoe and Mono Basin datasets
interpolated_x_tahoe = InterpolatedTahoeData.X;  % Replace 'x_values' with the actual field name in your data
interpolated_z_tahoe = InterpolatedTahoeData.Z;  % Replace 'z_values' with the actual field name in your data

fprintf('Analyzing Tahoe profile...\n');

%Store final profiles for Tahoe (trim the first 1550 points)
final_linear_diffusion_profile_tahoe = slope_height(1551:3000);  % From the linear diffusion model
final_grain_transport_profile_tahoe = Grain_transport_profile(1551:3000);  % From the grain transport model
final_non_local_profile_tahoe = z_nonlocal(1551:3000);  % From the non-local model
fprintf('Storing profiles\n');

%Load profiles
final_linear_diffusion_profile_tahoe = slope_height(1551:3000);  % From the linear diffusion model
final_grain_transport_profile_tahoe = Grain_transport_profile(1551:3000);  % From the grain transport model
final_non_local_profile_tahoe = z_nonlocal(1551:3000);  % From the non-local model

%Create an x-axis (assuming all profiles are of the same length)
x_values = 1551:(1550 + length(final_linear_diffusion_profile_tahoe));  % Adjust x-values to match original index range

%Interpolate Tahoe profiles to match the X values of the interpolated profile
x_model_tahoe = linspace(min(interpolated_x_tahoe), max(interpolated_x_tahoe), length(final_linear_diffusion_profile_tahoe));
linear_diffusion_profile_interp_tahoe = interp1(x_model_tahoe, final_linear_diffusion_profile_tahoe, interpolated_x_tahoe, 'linear');
grain_transport_profile_interp_tahoe = interp1(x_model_tahoe, final_grain_transport_profile_tahoe, interpolated_x_tahoe, 'linear');
non_local_profile_interp_tahoe = interp1(x_model_tahoe, final_non_local_profile_tahoe, interpolated_x_tahoe, 'linear');
fprintf('Interpolating profiles \n');

%Convert profiles to meters and normalize
linear_diffusion_profile_interp_tahoe = linear_diffusion_profile_interp_tahoe / 10;
grain_transport_profile_interp_tahoe = grain_transport_profile_interp_tahoe / 10;
non_local_profile_interp_tahoe = non_local_profile_interp_tahoe / 10;
interpolated_z_tahoe = interpolated_z_tahoe;
fprintf('Normalizing to meters\n');

linear_diffusion_profile_interp_tahoe = linear_diffusion_profile_interp_tahoe - min(linear_diffusion_profile_interp_tahoe);
grain_transport_profile_interp_tahoe = grain_transport_profile_interp_tahoe - min(grain_transport_profile_interp_tahoe);
non_local_profile_interp_tahoe = non_local_profile_interp_tahoe - min(non_local_profile_interp_tahoe);
interpolated_z_tahoe = interpolated_z_tahoe - min(interpolated_z_tahoe);
fprintf('Calculating differences\n');

%Calculate errors for Tahoe profile
rmse_grain_vs_tahoe = sqrt(mean((grain_transport_profile_interp_tahoe - interpolated_z_tahoe).^2));
rmse_nonlocal_vs_tahoe = sqrt(mean((non_local_profile_interp_tahoe - interpolated_z_tahoe).^2));
rmse_linear_diffusion_vs_tahoe = sqrt(mean((linear_diffusion_profile_interp_tahoe - interpolated_z_tahoe).^2));

mae_grain_vs_tahoe = mean(abs(grain_transport_profile_interp_tahoe - interpolated_z_tahoe));
mae_nonlocal_vs_tahoe = mean(abs(non_local_profile_interp_tahoe - interpolated_z_tahoe));
mae_linear_diffusion_vs_tahoe = mean(abs(linear_diffusion_profile_interp_tahoe - interpolated_z_tahoe));

max_abs_error_grain_vs_tahoe = max(abs(grain_transport_profile_interp_tahoe - interpolated_z_tahoe));
max_abs_error_nonlocal_vs_tahoe = max(abs(non_local_profile_interp_tahoe - interpolated_z_tahoe));
max_abs_error_linear_diffusion_vs_tahoe = max(abs(linear_diffusion_profile_interp_tahoe - interpolated_z_tahoe));


%Display the results for Tahoe profile
fprintf('Tahoe Profile Comparison:\n');
fprintf('Grain Transport Model - RMSE: %.4f, MAE: %.4f, Max Abs Error: %.4f\n', rmse_grain_vs_tahoe, mae_grain_vs_tahoe, max_abs_error_grain_vs_tahoe);
fprintf('Non-Local Model - RMSE: %.4f, MAE: %.4f, Max Abs Error: %.4f\n', rmse_nonlocal_vs_tahoe, mae_nonlocal_vs_tahoe, max_abs_error_nonlocal_vs_tahoe);
fprintf('Linear diffusion Model - RMSE: %.4f, MAE: %.4f, Max Abs Error: %.4f\n', rmse_linear_diffusion_vs_tahoe, mae_linear_diffusion_vs_tahoe, max_abs_error_linear_diffusion_vs_tahoe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           Plotting                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot linear diffusion, grain transport, and non-local models for
% comparison
figure(1) 
hold on
%axis equal
plot((1:6000) / 10, slope_original(1:6000) / 10, 'DisplayName', 'Initial Slope'); % Scale x and y by 10
plot((1:6000) / 10, slope_height(1:6000) / 10, 'DisplayName', 'Linear-diffusion model'); % Scale x and y by 10
plot((1:6000) / 10, Grain_transport_profile(1:6000) / 10, 'DisplayName', 'Grain-transport model'); % Scale x and y by 10
plot((1:6000) / 10, z_nonlocal(1:6000) / 10, 'DisplayName', 'Non-Local model');
ylim([0, 75])
title('Tahoe profile across full spatial domain');
xlabel('Horizontal Distance (m)'); % Horizontal distance in meters
ylabel('Height (m)'); % Vertical height in meters
legend();
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the interpolated profiles
    figure(2);
    hold on; % This allows multiple plots on the same figure

    plot(interpolated_x_tahoe, linear_diffusion_profile_interp_tahoe, 'DisplayName', 'Linear Diffusion');
    plot(interpolated_x_tahoe, grain_transport_profile_interp_tahoe, 'DisplayName', 'Grain Transport');
    plot(interpolated_x_tahoe, non_local_profile_interp_tahoe, 'DisplayName', 'Non-local Model');
    plot(interpolated_x_tahoe, interpolated_z_tahoe, 'DisplayName', 'Measured Hillslope Profile');

    xlabel('distance (m)');  % Label for x-axis
    ylabel('Elvation (m)');           % Label for y-axis
    title('Model profiles vs measured Tahoe profile');  % Title of the plot
    legend('location', 'northwest', 'FontSize', 11);  % Show the legend based on 'DisplayName'
    grid on;  % Add a grid for better readability

   % Adjust x-axis limits based on the actual range of interpolated_x_tahoe
        xlim([min(interpolated_x_tahoe), max(interpolated_x_tahoe)]);  % Auto-adjust x-axis to the interpolated range
    hold off;  % Release hold on the figure

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plotting the lowering of the crest over time for both models
figure(3);
hold on;
plot(time_nonlocal, crest_elevation_nonlocal / 10, 'DisplayName', 'Non-Local Model');
plot(time_diffusion, crest_elevation_diffusion / 10, 'DisplayName', 'Linear Diffusion Model');
plot(time_points, crest_height_over_time / 10, 'DisplayName', 'Grain Transport');
xlabel('Time (years)');
ylabel('Tahoe crest Elevation (m)');
ylim([40 80])
title('Lowering of the Tahoe Crest Over Time');
legend;
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allows profile to fit in distribution chart (no need to modify)
Grain_transport_profile = Grain_transport_profile / 10;

figure(4)
hold on
title('Tahoe moraine at 45 kyr', 'FontSize', 12)

%Set left y-axis for grain size distributions (Weight %)
yyaxis left

plot((1:length(active_percent_small)) / 10, active_percent_small, 'b-', 'LineStyle', '-', 'DisplayName', '< 1 mm'); % Scale x by 10
plot((1:length(active_percent_medium)) / 10, active_percent_medium, 'r-', 'LineStyle', '-', 'DisplayName', '1 - 2 mm'); % Scale x by 10
plot((1:length(active_percent_medlar)) / 10, active_percent_medlar, 'c-', 'LineStyle', '-', 'DisplayName', '2 - 8 mm'); % Scale x by 10
plot((1:length(active_percent_large)) / 10, active_percent_large, 'k-', 'LineStyle', '-', 'DisplayName', '8 - 19 mm'); % Scale x by 10
plot((1:length(active_percent_Very_Large)) / 10, active_percent_Very_Large, 'm-', 'LineStyle', '-', 'DisplayName', '>19 mm'); % Scale x by 10

%Set left y-axis color to black
ax = gca;
ax.YColor = 'k';

%Add a line representing the sum of all weight percentages
weight_sum = active_percent_small + active_percent_medium + ...
             active_percent_medlar + active_percent_large + active_percent_Very_Large;
plot((1:length(weight_sum)) / 10, weight_sum, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Total Weight %');

%Set labels and limits for left y-axis (Weight %)
xlabel('Horizontal Distance (m)', 'FontSize', 14) % Horizontal distance in meters
ylabel('Weight %', 'FontSize', 14)
ylim([-3 100])

%Add secondary vertical axis for elevation on the right
yyaxis right
plot((1:array_length) / 10, Grain_transport_profile(1:array_length) / 10, 'LineStyle', '--', 'Color',  [0.75 0.75 0.75], 'DisplayName', 'Grain-transport model elevation');
ylabel('Elevation (m)', 'FontSize', 14)

%Set y-axis color
ax.YColor = 'k';

%Adjust the limits for the right y-axis (Elevation)
ylim([min(Grain_transport_profile/10) max(Grain_transport_profile/10)])

%Legend and x-axis limits
legend('Location', 'northwest', 'FontSize', 11)
xlim([140 300]) % Scale x-limits by 10 (from decimeters to meters)

%Switch to right y-axis for elevation
yyaxis left

plot(290, 53, 'om', 'HandleVisibility', 'off');
plot(290, 1, 'ok', 'HandleVisibility', 'off');
plot(290, 3, 'oc', 'HandleVisibility', 'off');
plot(290, 5, 'or', 'HandleVisibility', 'off');
plot(290, 38, 'ob', 'HandleVisibility', 'off');


plot(226, 16, 'om', 'HandleVisibility', 'off');
plot(226, 8, 'ok', 'HandleVisibility', 'off');
plot(226, 15, 'oc', 'HandleVisibility', 'off');
plot(226, 24, 'or', 'HandleVisibility', 'off');
plot(226, 37, 'ob', 'HandleVisibility', 'off');

%157
plot(167, 1, 'om', 'HandleVisibility', 'off');
plot(167, 10, 'ok', 'HandleVisibility', 'off');
plot(167, 17, 'oc', 'HandleVisibility', 'off');
plot(167, 30, 'or', 'HandleVisibility', 'off');
plot(167, 42, 'ob', 'HandleVisibility', 'off');

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures for grain size distribution vs. time
x = [4500 4500];
y = [2 60];
x1 = [9990 9990];
y1 = [2 60];


figure(5)
hold on
title('Time vs. grain size distribution (footslope) for Tahoe', 'FontSize', 12)
set(gca, 'Color', 'w'); % Set axes background to white
set(gcf, 'Color', 'w'); % Set figure background to white
plot(small_footslope_array, 'b');
plot(medium_footslope_array, 'r');
plot(medlar_footslope_array, 'c');
plot(large_footslope_array, 'k');
plot(very_Large_footslope_array, 'm');
plot(x, y, 'LineStyle', '--', 'Color', [0.75 0.75 0.75]);
plot(x1, y1, 'LineStyle', '--', 'Color', [0.75 0.75 0.75]);

% Add a line representing the sum of all weight percentages
weight_sum_footslope = small_footslope_array + medium_footslope_array + medlar_footslope_array + large_footslope_array + very_Large_footslope_array;
plot(weight_sum_footslope, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Total Weight %');

xticks([10 2000 4000 6000 8000 10000]);
xticklabels({'0', '20', '40', '60', '80', '100'});
xlabel('Time (kyr)', 'FontSize', 14);
ylabel('Weight %', 'FontSize', 14);
legend({'< 1 mm', '1 - 2 mm', '2 - 8 mm', '8 - 19 mm', '> 19 mm', 'Total Weight %'}, 'Location', 'northwest', 'FontSize', 11);
xlim([1 10000])
ylim([-3 100])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
hold on
title('Time vs. grain size distribution (Midslope) for Tahoe', 'FontSize', 12)
plot(small_midslope_array, 'b');
plot(medium_midslope_array, 'r');
plot(medlar_midslope_array, 'c');
plot(large_midslope_array, 'k');
plot(very_Large_midslope_array, 'm');
plot(x, y, 'LineStyle', '--', 'Color', [0.75 0.75 0.75]);
plot(x1, y1, 'LineStyle', '--', 'Color', [0.75 0.75 0.75]);

% Add a line representing the sum of all weight percentages
weight_sum_crest = small_crest_array + medium_crest_array + ...
                   medlar_crest_array + large_crest_array + very_Large_crest_array;
plot(weight_sum_crest, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Total Weight %');

xticks([10 2000 4000 6000 8000 10000]);
xticklabels({'0', '20', '40', '60', '80', '100'});
xlabel('Time (kyr)', 'FontSize', 14);
ylabel('Weight %', 'FontSize', 14);
legend({'< 1 mm', '1 - 2 mm', '2 - 8 mm', '8 - 19 mm', '> 19 mm', 'Total Weight %'}, 'Location', 'northwest', 'FontSize', 11);
xlim([1 10000])
ylim([-3 100])
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7)
hold on
title('Time vs. grain size distribution (crest) for Tahoe', 'FontSize', 12)
plot(small_crest_array, 'b');
plot(medium_crest_array, 'r');
plot(medlar_crest_array, 'c');
plot(large_crest_array, 'k');
plot(very_Large_crest_array, 'm');
plot(x, y, 'LineStyle', '--', 'Color', [0.75 0.75 0.75]);
plot(x1, y1, 'LineStyle', '--', 'Color', [0.75 0.75 0.75]);

% Add a line representing the sum of all weight percentages
weight_sum_crest = small_crest_array + medium_crest_array + ...
                   medlar_crest_array + large_crest_array + very_Large_crest_array;
plot(weight_sum_crest, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Total Weight %');

xticks([10 2000 4000 6000 8000 10000]);
xticklabels({'0', '20', '40', '60', '80', '100'});
xlabel('Time (kyr)', 'FontSize', 14);
ylabel('Weight %', 'FontSize', 14);
legend({'< 1 mm', '1 - 2 mm', '2 - 8 mm', '8 - 19 mm', '> 19 mm', 'Total Weight %'}, 'Location', 'northwest', 'FontSize', 11);
xlim([1 10000])
ylim([-3 100])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8);
% 1/200 since the model was originally written in decimeters and only one half of the hillslope is calculated.
small_footslope_rate = small_footslope_rate / 200;
small_midslope_rate = small_midslope_rate / 200;
small_crest_rate = small_crest_rate / 200;
plot(time_array_rates, small_footslope_rate, '-r', 'DisplayName', 'Footslope');
hold on;
plot(time_array, small_midslope_rate, '-g', 'DisplayName', 'Midslope');
plot(time_array, small_crest_rate, '-b', 'DisplayName', 'Crest');
xlabel('Time (years)');
ylabel('Transport Rate (qSmall)');
title('Transport Rate for <1 mm Grain Size');
legend;
hold off;

figure(9);
% Halving the rates since the true diffusion constant must be halved.
medium_footslope_rate = medium_footslope_rate / 200;
medium_midslope_rate = medium_midslope_rate / 200;
medium_crest_rate = medium_crest_rate / 200;
plot(time_array_rates, medium_footslope_rate, '-r', 'DisplayName', 'Footslope');
hold on;
plot(time_array_rates, medium_midslope_rate, '-g', 'DisplayName', 'Midslope');
plot(time_array_rates, medium_crest_rate, '-b', 'DisplayName', 'Crest');
xlabel('Time (years)');
ylabel('Transport Rate (qMedium)');
title('Transport Rate for 1-2 mm Grain Size');
legend;
hold off;

figure(10);
% Halving the rates since the true diffusion constant must be halved.
medlar_footslope_rate = medlar_footslope_rate / 200;
medlar_midslope_rate = medlar_midslope_rate / 200;
medlar_crest_rate = medlar_crest_rate / 200;
plot(time_array_rates, medlar_footslope_rate, '-r', 'DisplayName', 'Footslope');
hold on;
plot(time_array_rates, medlar_midslope_rate, '-g', 'DisplayName', 'Midslope');
plot(time_array_rates, medlar_crest_rate, '-b', 'DisplayName', 'Crest');
xlabel('Time (years)');
ylabel('Transport Rate (qMedlar)');
title('Transport Rate for 2-8 mm Grain Size');
legend;
hold off;


figure(11);
% Halving the rates since the true diffusion constant must be halved.
large_footslope_rate = large_footslope_rate / 200;
large_midslope_rate = large_midslope_rate / 200;
large_crest_rate = large_crest_rate / 200;
plot(time_array_rates, large_footslope_rate, '-r', 'DisplayName', 'Footslope');
hold on;
plot(time_array_rates, large_midslope_rate, '-g', 'DisplayName', 'Midslope');
plot(time_array_rates, large_crest_rate, '-b', 'DisplayName', 'Crest');
xlabel('Time (years)');
ylabel('Transport Rate (qLarge)');
title('Transport Rate for 8-19mm Grain Size');
legend;
hold off;

figure(12);
% Halving the rates since the true diffusion constant must be halved.
veryLarge_footslope_rate = veryLarge_footslope_rate / 200;
veryLarge_midslope_rate = veryLarge_midslope_rate / 200;
veryLarge_crest_rate = veryLarge_crest_rate / 200;
plot(time_array_rates, veryLarge_footslope_rate, '-r', 'DisplayName', 'Footslope');
hold on;
plot(time_array_rates, veryLarge_midslope_rate, '-g', 'DisplayName', 'Midslope');
plot(time_array_rates, veryLarge_crest_rate, '-b', 'DisplayName', 'Crest');
xlabel('Time (years)');
ylabel('Transport Rate (qVLarge)');
title('Transport Rate for >19 mm Grain Size');
legend;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Code Functions                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: Calculate RMSE for a given set of expected and observed values
function RMSE = calculate_rmse(expected_values, observed_values)
    num_points = length(expected_values); % Calculate number of data points
    squared_errors = 0; % Initialize sum of squared errors
    for i = 1:num_points
        squared_errors = squared_errors + (expected_values(i) - observed_values(i))^2;
    end
    RMSE = sqrt(squared_errors / num_points); % Calculate RMSE
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: Calculate volume flux for a given grain size
function q = calculate_flux(K, dydx, dt, dx)
    q = K * dydx * (dt / dx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: Calculate grain movement based on slope angles
function [added_volume, lost_volume] = calculate_grain_movement(slope_angle, qx, added_volume, lost_volume)
    array_length = length(slope_angle); % Get the length of the slope_angle array
    for j = array_length:-1:5
        if slope_angle(j) >= 0
            %Add volume to the downslope segment and subtract from the current one
            added_volume(j - 1) = added_volume(j - 1) + qx(j);
            lost_volume(j) = -qx(j); % Lost volume is negative of the flux
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function: Calculate fluxes for each grain size class

function [qxSmall, qxMedium, qxMedLar, qxLarge, qxVeryLarge, qxTotal] = ...
    calculate_fluxes(array_length, qTotal,  qSmall, qMedium, qMedLar, qLarge, qVeryLarge, ...
                     active_percent_small, active_percent_medium, ...
                     active_percent_medlar, active_percent_large, active_percent_Very_Large)

    %Initialize output variables
    qxSmall = zeros(1, array_length);
    qxMedium = zeros(1, array_length);
    qxMedLar = zeros(1, array_length);
    qxLarge = zeros(1, array_length);
    qxVeryLarge = zeros(1, array_length);
    qxTotal = zeros(1, array_length);

    for i = 5:array_length
        if qTotal(i) > 0
            %Calculate qz for each grain size
            qzSmall = qSmall(i) * active_percent_small(i);
            qzMedium = qMedium(i) * active_percent_medium(i);
            qzMedLar = qMedLar(i) * active_percent_medlar(i);
            qzLarge = qLarge(i) * active_percent_large(i);
            qzVeryLarge = qVeryLarge(i) * active_percent_Very_Large(i);
            qzTotal = qzSmall + qzMedium + qzMedLar + qzLarge + qzVeryLarge;

            if qzTotal > 0
                %Recalculate qx for each grain size based on total
                qxSmall(i) = (qzSmall / qzTotal) * qTotal(i);
                qxMedium(i) = (qzMedium / qzTotal) * qTotal(i);
                qxMedLar(i) = (qzMedLar / qzTotal) * qTotal(i);
                qxLarge(i) = (qzLarge / qzTotal) * qTotal(i);
                qxVeryLarge(i) = (qzVeryLarge / qzTotal) * qTotal(i);
                qxTotal(i) = qxSmall(i) + qxMedium(i) + qxMedLar(i) + qxLarge(i) + qxVeryLarge(i);
            else
                %If qzTotal is zero or negative, skip recalculation
                qxSmall(i) = 0;
                qxMedium(i) = 0;
                qxMedLar(i) = 0;
                qxLarge(i) = 0;
                qxVeryLarge(i) = 0;
                qxTotal(i) = 0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Optimized sediment flux function for the non-local model
function q_non_local = sediment_flux(z, weights, dx, K_nonlocal)
    %Local gradient of the slope
    gradient_local = gradient(z, dx);

    %Calculate non-local influence using matrix multiplication with precomputed weights
    grad_z = gradient(z, dx);  % Local gradient vector
    non_local_gradient = weights * grad_z';  % Sparse matrix multiplication for non-local influence

    %Total flux is the sum of local and non-local components
    q_non_local = -K_nonlocal * (gradient_local + non_local_gradient');  % Assign total flux to q_non_local
end