%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      LHS Model                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tStart = tic;

% Duration of simulation
maxtime = 45000; % Years (Tahoe: 45000; Mono Basin: 100000)
array_length = 6000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Grain transport & distribution diffusion model           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time increment
dt = 0.2; % Years

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

% Initial grain percents- Taken from Mono Basin crest. The assumption is
% that there is little/no transport directly at the crest.
initial_percent_small = 0.42
initial_percent_medium = 0.12; 
initial_percent_medlar = 0.12; 
initial_percent_large = 0.14; 
initial_percent_Very_Large = 0.2; 


% Define parameter ranges for LHS
param_ranges = [
    0.0005, 0.04;    % Range for KSmall
    0.0005, 0.04;    % Range for KMedium
    0.0005, 0.04;    % Range for KMedLar
    0.0005, 0.04;    % Range for KLarge
    0.0005, 0.04;    % Range for KVeryLarge
    0.000001, 0.00004; % Range for very_Large_wr
    0.0000001, 0.000005; % Range for large_wr
    0.0000001, 0.000005; % Range for medlar_wr
    0.0000001, 0.000005; % Range for medium_wr
];

target_sum = 0.031; % Target sum for the k variables

% Expected values for Tahoe moraine (XSmall, Small, Medium, MedLar, Large)
% [footslope, midslope, crest]
expected_small = [24, 32, 47];
expected_medium = [17, 25, 7];
expected_medlar = [37, 22, 6];
expected_large = [20, 8, 17];
expected_Very_Large = [2, 13, 23];

num_samples = 4; % Set the number of LHS samples
valid_samples = []; % Array to store valid samples

% Loop until we have the required number of valid samples
while size(valid_samples, 1) < num_samples
    % Generate a single LHS sample
    lhs_samples = lhsdesign(1, size(param_ranges, 1)); % Generate one sample at a time
    param_set = lhs_samples .* (param_ranges(:,2)' - param_ranges(:,1)') + param_ranges(:,1)';

    % Scale the k values to ensure their sum is equal to target_sum
    k_values = param_set(1:5);
    k_total = sum(k_values);
    scaled_k_values = (k_values / k_total) * target_sum;
    param_set(1:5) = scaled_k_values; % Replace original k values with scaled ones

    % Extract weathering parameters
    very_Large_wr = param_set(6);
    large_wr = param_set(7);
    medlar_wr = param_set(8);
    medium_wr = param_set(9);

    % Check if the constraints are satisfied
    if very_Large_wr > large_wr && large_wr > medlar_wr && medlar_wr > medium_wr
        % Add small_wr as zero and append to valid_samples
        param_set = [param_set, 0]; % Ensure 10th parameter (small_wr) is included
        valid_samples = [valid_samples; param_set];
    end
end

% Initialize arrays to store RMSE and KTotal values
rmse_values = zeros(num_samples, 1);
k_total_values = zeros(num_samples, 1);

% Calculate the number of time steps when rem(time,10) == 0
num_appends = floor(maxtime / 10) + 1; % Adding 1 to include the initial time (time = 0)



 %Main LHS loop - can be parallelized using parfor
parfor sample_idx = 1:num_samples

    % Fetch the valid sample for the current worker
    params = valid_samples(sample_idx, :);


    % Initialize the grain sizes for this iteration to prevent reduction variable error
    initial_percent_small = 0.38; 
    initial_percent_medium = 0.06; 
    initial_percent_medlar = 0.01; 
    initial_percent_large = 0.05; 
    initial_percent_Very_Large = 0.50;

    KSmall = params(1);
    KMedium = params(2);
    KMedLar = params(3);
    KLarge = params(4);
    KVeryLarge = params(5);
    very_Large_wr = params(6);
    large_wr = params(7);
    medlar_wr = params(8);
    medium_wr = params(9);
    small_wr = params(10);

    % Compute KTotal and store it
    KTotal = sum(params(1:5));
    k_total_values(sample_idx) = KTotal;

    % Adjust initial grain sizes for each sample
    adjusted_percent_small = initial_percent_small * (1 + medium_wr);
    adjusted_percent_medium = initial_percent_medium * (1 - medium_wr);
    adjusted_percent_medlar = initial_percent_medlar * (1 - medlar_wr);
    adjusted_percent_large = initial_percent_large * (1 - large_wr);
    adjusted_percent_Very_Large = initial_percent_Very_Large * (1 - very_Large_wr);

    % Assign adjusted values to active percentages
    active_percent_small = ones(1, array_length) * adjusted_percent_small;
    active_percent_medium = ones(1, array_length) * adjusted_percent_medium;
    active_percent_medlar = ones(1, array_length) * adjusted_percent_medlar;
    active_percent_large = ones(1, array_length) * adjusted_percent_large;
    active_percent_Very_Large = ones(1, array_length) * adjusted_percent_Very_Large;
    

    % Initialize profile and perform main simulation loop
    xdistance = zeros(1, array_length);
    ydistance = zeros(1, array_length);
    alfa = 28;

    for counter = 1:1800
        xdistance(counter) = counter;
        ydistance(counter) = 0;
    end
    for counter = 1801:3000
        xdistance(counter) = counter;
        ydistance(counter) = (counter - 1800) * tan((alfa / 360) * 2 * pi);
    end
    for counter = 3001:4200
        xdistance(counter) = counter;
        ydistance(counter) = ydistance(6001 - counter);
    end
    for counter = 4201:array_length
        xdistance(counter) = counter;
        ydistance(counter) = 0;
    end

    slope_original = ydistance;
    node_count = length(xdistance);

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
    
    save_interval = 10; % Save the profile every 10 years
    save_idx = 1; % Index for storing profiles
    
    % Initialize index for appending data
    idx = 1;
    
    % Imports raw profile to now use
    Grain_transport_profile = slope_original;
    
    dydx = zeros(1, array_length);
    for i = 2:(array_length - 1)
        dydx(i) = (Grain_transport_profile(i + 1) - Grain_transport_profile(i - 1)) / (2 * dx);
    end
    
    fprintf('Begin Grain Transport and distribution Simulation\n');
    
    %%%%%%%%%%%%%%%%%%%%
    
    % Main simulation loop
    for time = 0:dt:maxtime
    
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
                % Calculate new volumes with safeguards for each grain size
                new_volume_small(i) = ((active_volume - net_volume_total(i)) * active_percent_small(i)) + net_volume_small(i);
                new_volume_medium(i) = ((active_volume - net_volume_total(i)) * active_percent_medium(i)) + net_volume_medium(i);
                new_volume_medlar(i) = ((active_volume - net_volume_total(i)) * active_percent_medlar(i)) + net_volume_medlar(i);
                new_volume_large(i) = ((active_volume - net_volume_total(i)) * active_percent_large(i)) + net_volume_large(i);
                new_volume_Very_Large(i) = ((active_volume - net_volume_total(i)) * active_percent_Very_Large(i)) + net_volume_Very_Large(i);
            
                % Calculate the total volume and normalize
                new_volume_total(i) = new_volume_small(i) + new_volume_medium(i) + new_volume_medlar(i) + new_volume_large(i) + new_volume_Very_Large(i);
            
                % If total volume is valid, proceed with normalization
                if new_volume_total(i) > 0 && isfinite(new_volume_total(i))
                    active_percent_small(i) = max(0, min(1, new_volume_small(i) / new_volume_total(i)));
                    active_percent_medium(i) = max(0, min(1, new_volume_medium(i) / new_volume_total(i)));
                    active_percent_medlar(i) = max(0, min(1, new_volume_medlar(i) / new_volume_total(i)));
                    active_percent_large(i) = max(0, min(1, new_volume_large(i) / new_volume_total(i)));
                    active_percent_Very_Large(i) = max(0, min(1, new_volume_Very_Large(i) / new_volume_total(i)));
            
                    % Additional normalization
                    total_percent = active_percent_small(i) + active_percent_medium(i) + active_percent_medlar(i) + active_percent_large(i) + active_percent_Very_Large(i);
                    if abs(total_percent - 1) > 1e-6
                        active_percent_small(i) = active_percent_small(i) / total_percent;
                        active_percent_medium(i) = active_percent_medium(i) / total_percent;
                        active_percent_medlar(i) = active_percent_medlar(i) / total_percent;
                        active_percent_large(i) = active_percent_large(i) / total_percent;
                        active_percent_Very_Large(i) = active_percent_Very_Large(i) / total_percent;
                    end
                else
                    % Debug output if new_volume_total is invalid
                    fprintf('Invalid new_volume_total at index %d: %f\n', i, new_volume_total(i));
                end
            
                % Debugging output to track active_percent_* values
                if any([active_percent_small(i), active_percent_medium(i), active_percent_medlar(i), active_percent_large(i), active_percent_Very_Large(i)] < 0) || ...
                   any([active_percent_small(i), active_percent_medium(i), active_percent_medlar(i), active_percent_large(i), active_percent_Very_Large(i)] > 1)
                    fprintf('Active_percent values out of range at index %d\n', i);
                    fprintf('Small: %f, Medium: %f, MedLar: %f, Large: %f, Very_Large: %f\n', ...
                            active_percent_small(i), active_percent_medium(i), active_percent_medlar(i), ...
                            active_percent_large(i), active_percent_Very_Large(i));
                end
            
                % Final assertion
                assert(abs(active_percent_small(i) + active_percent_medium(i) + active_percent_medlar(i) + active_percent_large(i) + active_percent_Very_Large(i) - 1) < 1e-6, ...
                       'Error: The sum of the percentages in accumulation does not equal 1 at index %d', i);
            end


    
             % Weathering
             new_w_volume_Very_Large(i) = new_volume_Very_Large(i) * (1 - very_Large_wr);
             new_w_volume_large(i) = new_volume_large(i) + (very_Large_wr * new_volume_Very_Large(i));
             new_w_volume_medlar(i) = new_volume_medlar(i) + (large_wr * new_volume_large(i));
             new_w_volume_medium(i) = new_volume_medium(i) + (medlar_wr * new_w_volume_medlar(i));
             new_w_volume_small(i) = new_volume_small(i) + (medium_wr * new_w_volume_medium(i));
             new_w_volume_total(i) = new_w_volume_small(i) + new_w_volume_medium(i) + new_w_volume_medlar(i) + new_w_volume_large(i) + new_w_volume_Very_Large(i);
                
             if new_w_volume_total(i) > 0
                 active_percent_small(i) = new_w_volume_small(i) / new_w_volume_total(i);
                 active_percent_medium(i) = new_w_volume_medium(i) / new_w_volume_total(i);
                 active_percent_medlar(i) = new_w_volume_medlar(i) / new_w_volume_total(i);
                 active_percent_large(i) = new_w_volume_large(i) / new_w_volume_total(i);
                 active_percent_Very_Large(i) = new_w_volume_Very_Large(i) / new_w_volume_total(i);
             end
                
             % Additional normalization to correct any minor floating-point deviations
             total_percent = active_percent_small(i) + active_percent_medium(i) + ...
                                active_percent_medlar(i) + active_percent_large(i) + active_percent_Very_Large(i);
                
             % Adjust values if the sum deviates from 1
             if abs(total_percent - 1) > 1e-6  % Allow small tolerance for floating-point precision
                  active_percent_small(i) = active_percent_small(i) / total_percent;
                  active_percent_medium(i) = active_percent_medium(i) / total_percent;
                  active_percent_medlar(i) = active_percent_medlar(i) / total_percent;
                  active_percent_large(i) = active_percent_large(i) / total_percent;
                  active_percent_Very_Large(i) = active_percent_Very_Large(i) / total_percent;
             end
                
                % Debug output for inspection if an error is triggered
             if abs(active_percent_small(i) + active_percent_medium(i) + active_percent_medlar(i) + ...
                       active_percent_large(i) + active_percent_Very_Large(i) - 1) > 1e-6
                    fprintf('Normalization error at index %d\n', i);
                    fprintf('Sum of percentages: %f\n', total_percent);
                    fprintf('Small: %f, Medium: %f, MedLar: %f, Large: %f, Very_Large: %f\n', ...
                            active_percent_small(i), active_percent_medium(i), active_percent_medlar(i), ...
                            active_percent_large(i), active_percent_Very_Large(i));
         end
                
                % Final assertion to catch any remaining discrepancies
                assert(abs(active_percent_small(i) + active_percent_medium(i) + ...
                           active_percent_medlar(i) + active_percent_large(i) + ...
                           active_percent_Very_Large(i) - 1) < 1e-6, ...
                       'Error: The sum of the percentages in weathering does not equal 1 at index %d', i)


    

    
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
            small_footslope = active_percent_small(1570) * 100; % Footslope = 1860 for Tahoe
            medium_footslope = active_percent_medium(1570) * 100;
            medlar_footslope = active_percent_medlar(1570) * 100;
            large_footslope = active_percent_large(1570) * 100;
            very_Large_footslope = active_percent_Very_Large(1570) * 100;
    
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
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%          Begin Grain Size Distribution Error Calculations           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Convert active percent to out of 100 for aesthetics
    active_percent_small = active_percent_small * 100;
    active_percent_medium = active_percent_medium * 100;
    active_percent_medlar = active_percent_medlar * 100;
    active_percent_large = active_percent_large * 100;
    active_percent_Very_Large = active_percent_Very_Large * 100;
    
    % Main script for RMSE calculations
    % Tahoe moraine data points (Foot, Mid, Crest) for each grain size
    points = [1570, 2400, 2940]; % Location based on % distance from crest whereby the bulk sample grain size distributions were calculated
    
    % Observed values for Tahoe moraine 
    observed_small = [active_percent_small(1570), active_percent_small(2400), active_percent_small(2940)];
    observed_medium = [active_percent_medium(1570), active_percent_medium(2400), active_percent_medium(2940)];
    observed_medlar = [active_percent_medlar(1570), active_percent_medlar(2400), active_percent_medlar(2940)];
    observed_large = [active_percent_large(1570), active_percent_large(2400), active_percent_large(2940)];
    observed_Very_Large = [active_percent_Very_Large(1570), active_percent_Very_Large(2400), active_percent_Very_Large(2940)];


    % Calculate RMSE for each grain size
    RMSE_Small_Young = calculate_rmse(expected_small, observed_small);
    RMSE_Medium_Young = calculate_rmse(expected_medium, observed_medium);
    RMSE_MedLar_Young = calculate_rmse(expected_medlar, observed_medlar);
    RMSE_Large_Young = calculate_rmse(expected_large, observed_large);
    RMSE_Very_Large_Young = calculate_rmse(expected_Very_Large, observed_Very_Large);

    % Total RMSE for Tahoe moraine
    RMSE_Total_Young = (RMSE_Small_Young + RMSE_Medium_Young + ...
                            RMSE_MedLar_Young + RMSE_Large_Young + RMSE_Very_Large_Young) / 5;

    % Store the RMSE value for this parameter set
    rmse_values(sample_idx) = RMSE_Total_Young;  

end

% Write results to CSV
output_filename = 'parameter_rmse_data.csv';
data_table = array2table([valid_samples, k_total_values, rmse_values], ...
    'VariableNames', {'KSmall', 'KMedium', 'KMedLar', 'KLarge', 'KVeryLarge', ...
                      'very_Large_wr', 'large_wr', 'medlar_wr', 'medium_wr', 'small_wr', ...
                      'KTotal', 'RMSE_Total_Young'});
if isfile(output_filename)
    writetable(data_table, output_filename, 'WriteMode', 'append');
else
    writetable(data_table, output_filename);
end

tEnd = toc(tStart);
fprintf('Simulation completed in %.2f seconds.\n', tEnd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Code Functions                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: Calculate RMSE for a given set of expected and observed values
function RMSE = calculate_rmse(expected_values, observed_values)
    num_points = length(expected_values); % Calculate number of data points
    squared_errors = 0; % Initialize sum of squared errors
    for i = 1:num_points
        squared_errors = squared_errors + (expected_values(i) - observed_values(i))^2;
    end
    RMSE = sqrt(squared_errors / num_points); % Calculate RMSE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: Calculate volume flux for a given grain size
function q = calculate_flux(K, dydx, dt, dx)
    q = K * dydx * (dt / dx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: Calculate grain movement based on slope angles
function [added_volume, lost_volume] = calculate_grain_movement(slope_angle, qx, added_volume, lost_volume)
    array_length = length(slope_angle); % Get the length of the slope_angle array
    for j = array_length:-1:5
        if slope_angle(j) >= 0
            % Add volume to the downslope segment and subtract from the current one
            added_volume(j - 1) = added_volume(j - 1) + qx(j);
            lost_volume(j) = -qx(j); % Lost volume is negative of the flux
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: Calculate fluxes for each grain size class
function [qxSmall, qxMedium, qxMedLar, qxLarge, qxVeryLarge, qxTotal] = ...
    calculate_fluxes(array_length, qTotal,  qSmall, qMedium, qMedLar, qLarge, qVeryLarge, ...
                     active_percent_small, active_percent_medium, ...
                     active_percent_medlar, active_percent_large, active_percent_Very_Large)

    % Initialize output variables
    qxSmall = zeros(1, array_length);
    qxMedium = zeros(1, array_length);
    qxMedLar = zeros(1, array_length);
    qxLarge = zeros(1, array_length);
    qxVeryLarge = zeros(1, array_length);
    qxTotal = zeros(1, array_length);

    for i = 5:array_length
        if qTotal(i) > 0
            % Calculate qz for each grain size
            qzSmall = qSmall(i) * active_percent_small(i);
            qzMedium = qMedium(i) * active_percent_medium(i);
            qzMedLar = qMedLar(i) * active_percent_medlar(i);
            qzLarge = qLarge(i) * active_percent_large(i);
            qzVeryLarge = qVeryLarge(i) * active_percent_Very_Large(i);
            qzTotal = qzSmall + qzMedium + qzMedLar + qzLarge + qzVeryLarge;

            if qzTotal > 0
                % Recalculate qx for each grain size based on total
                qxSmall(i) = (qzSmall / qzTotal) * qTotal(i);
                qxMedium(i) = (qzMedium / qzTotal) * qTotal(i);
                qxMedLar(i) = (qzMedLar / qzTotal) * qTotal(i);
                qxLarge(i) = (qzLarge / qzTotal) * qTotal(i);
                qxVeryLarge(i) = (qzVeryLarge / qzTotal) * qTotal(i);
                qxTotal(i) = qxSmall(i) + qxMedium(i) + qxMedLar(i) + qxLarge(i) + qxVeryLarge(i);
            else
                % If qzTotal is zero or negative, skip recalculation
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
