function [source_position, k_opt, theta_opt, phi_opt] = estimate_voltage_source_gd(positions, voltages)
    % Estimate dipole source location and parameters using gradient descent
    
    % Parameters
    alpha = 0.005;  % Adaptive learning rate
    beta = 0.9;  % Momentum factor
    max_iters = 5000; % Max iterations
    tol = 1e-8;  % Convergence tolerance
    epsilon = 1e-6; % Small step for numerical gradient

    % Normalize voltages for stable loss function
    voltage_magnitude = max(abs(voltages));
    normalized_voltages = voltages / voltage_magnitude;

    % Multiple random initializations to escape local minima
    best_loss = inf;
    best_vars = [];

    for restart = 1:3  % Try 3 different initializations
        % Randomized initial guess within a reasonable range
        vars = [mean(positions) + randn(1,3) * 2, 1 + randn, pi/4 + randn * 0.1, pi/4 + randn * 0.1];

        % Initialize momentum
        velocity = zeros(size(vars));

        % Run gradient descent
        for iter = 1:max_iters
            prev_vars = vars; % Store previous values
            
            % Compute gradients for all parameters at once
            grads = zeros(size(vars));
            for i = 1:length(vars)
                grads(i) = compute_gradient(positions, normalized_voltages, vars, i, epsilon);
            end
            
            % Update parameters with momentum
            velocity = beta * velocity - alpha * grads;  
            vars = vars + velocity;

            % Compute loss for this iteration
            current_loss = compute_loss(positions, normalized_voltages, vars);

            % Keep track of the best solution found
            if current_loss < best_loss
                best_loss = current_loss;
                best_vars = vars;
            end

            % Check convergence
            if norm(vars - prev_vars) < tol
                break;
            end
        end
    end

    % Extract optimized parameters
    source_position = best_vars(1:3);  % [x_s, y_s, z_s]
    k_opt = best_vars(4);  % Dipole strength constant
    theta_opt = best_vars(5);  % Polar angle
    phi_opt = best_vars(6);  % Azimuthal angle
end

function loss = compute_loss(positions, voltages, vars)
    % Computes the relative error loss
    predicted_voltages = model_dipole_voltage(positions, vars(1:3), vars(4), vars(5), vars(6));
    loss = sum(((voltages - predicted_voltages) ./ voltages).^2);
end

function grad = compute_gradient(positions, voltages, vars, param_idx, epsilon)
    % Compute numerical gradient using central difference method
    vars_plus = vars;
    vars_minus = vars;
    vars_plus(param_idx) = vars_plus(param_idx) + epsilon;
    vars_minus(param_idx) = vars_minus(param_idx) - epsilon;

    % Compute function values at perturbed points
    loss_plus = compute_loss(positions, voltages, vars_plus);
    loss_minus = compute_loss(positions, voltages, vars_minus);

    % Compute gradient approximation
    grad = (loss_plus - loss_minus) / (2 * epsilon);
end

function modeled_voltages = model_dipole_voltage(positions, source, scalar_constant, theta, phi)
    % Compute modeled voltage of a dipole
    r_vectors = positions - source;  % Nx3 matrix of position vectors from source to points
    distances = sqrt(sum(r_vectors.^2, 2));  % Euclidean distances

    % Unit vectors from source to measurement points
    unit_r_vectors = r_vectors ./ distances;

    % Dipole direction vector
    dipole_direction = [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)];

    % Compute dot product between unit_r_vectors and dipole direction
    dot_products = sum(unit_r_vectors .* dipole_direction, 2);

    % Compute modeled voltage using inverse-square law
    modeled_voltages = scalar_constant .* dot_products ./ (distances.^2);
end
