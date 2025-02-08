%% optimize the source location and the constant multiple of the field to use

function [source_position, k_opt, theta_opt, phi_opt, modeled_voltages] = estimate_voltage_source_v3(positions, voltages)
    % Function to estimate the source of a dipole's voltage in space along with constant k, theta, and phi
    % Input:
    % - positions: Nx3 matrix where each row represents the [x, y, z] coordinates of the measurement point
    % - voltages: Nx1 vector where each element represents the voltage at the corresponding measurement point
    % Output:
    % - source_position: Estimated [x, y, z] coordinates of the dipole source
    % - k_opt: Optimized constant k (related to the dipole moment)
    % - theta_opt: Optimized dipole angle theta (polar angle with respect to z-axis)
    % - phi_opt: Optimized dipole angle phi (azimuthal angle in the xy-plane)
    % - modeled_voltages: The computed modeled voltages using the optimized parameters

    % Number of measurement points
    N = size(positions, 1);

    % Normalize the input voltages for relative error calculation
    voltage_magnitude = max(abs(voltages));
    normalized_voltages = voltages / voltage_magnitude;

    % Initial guess for the source position [x_s, y_s, z_s], k, theta, and phi
    initial_guess = [mean(positions), 1, pi/4, pi/4];  % [x_s, y_s, z_s, k, theta, phi]

    % Objective function to minimize: relative error between measured and modeled voltages
    obj_fun = @(vars) sum(((normalized_voltages - model_dipole_voltage(positions, vars(1:3), vars(4), vars(5), vars(6))) ./ normalized_voltages).^2);

    % Use fminsearch to find the source position, k, theta, and phi that minimize the objective function
    options = optimset('Display', 'off');  % Turn off display for fminsearch
    optimal_vars = fminsearch(obj_fun, initial_guess, options);

    % Extract the optimized source position, k value, theta, and phi
    source_position = optimal_vars(1:3);  % [x_s, y_s, z_s]
    k_opt = optimal_vars(4);  % optimized k (constant)
    theta_opt = optimal_vars(5);  % optimized theta (polar angle)
    phi_opt = optimal_vars(6);  % optimized phi (azimuthal angle)

    % Calculate the modeled voltages using the optimized parameters
    modeled_voltages = model_dipole_voltage(positions, source_position, k_opt, theta_opt, phi_opt) * voltage_magnitude; % Scale back to original units
end

function modeled_voltages = model_dipole_voltage(positions, source_position, scalar_constant, theta, phi)
    % Function to compute modeled voltage of a dipole based on its orientation and position
    % positions: Nx3 matrix of measurement positions
    % source: 1x3 vector representing the [x, y, z] coordinates of the source
    % k: constant related to dipole moment
    % theta: polar angle of the dipole (with respect to the z-axis)
    % phi: azimuthal angle of the dipole (in the xy-plane)
    % Returns: modeled_voltages, an Nx1 vector of computed voltages

    % Compute distances between each measurement point and the source
    r_vectors = positions - source_position;  % Nx3 matrix of position vectors from source to points
    distances = sqrt(sum(r_vectors.^2, 2));  % Euclidean distances from source to each point

    % Unit vectors pointing from the source to each measurement point
    unit_r_vectors = r_vectors ./ distances;

    % Dipole direction unit vector, parameterized by theta and phi
    dipole_direction = [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)];

    % Compute the dot product between unit_r_vectors and the dipole direction
    dot_products = sum(unit_r_vectors .* dipole_direction, 2);

    % Calculate the modeled voltage using the dipole potential formula (inverse-square law with angular dependency)
    modeled_voltages = scalar_constant .* dot_products ./ (distances.^2);
end
