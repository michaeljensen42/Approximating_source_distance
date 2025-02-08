function [time_series_all, srate, t_tot, interval] = generate_BB_time_series_v4(frate_high, frate_low, trials, units)

% Input parameters
srate = 10000; % Sampling rate (10 kHz)
lat_s = 0.002; % Latency before decay kicks in (seconds)
t_tot = (trials * 3) + 3; % Total duration (seconds)
interval = 3; % Duration for each epoch (3 seconds)
psn_num = 6000; % Number of presynaptic neurons (6000 for a pyramidal neuron)
t = 0:floor(0.05 * srate); % 50 ms window for post-synaptic potential

epochs = t_tot / interval; % Total number of epochs

high_duration = 0.1; % Duration of high firing rate in seconds
high_samples = floor(high_duration * srate); % Number of steps for high firing rate

for un = 1:units

    time_series = zeros(1, t_tot * srate); % Preallocate time series

    % Loop through each epoch
    for epoch = 2:epochs
    
        % Determine if this epoch has a high firing rate for the first 200 ms
        is_high = 1; % Stimulate every epoch
        
        t_epoch = interval; % Time for this epoch (3 seconds)
        dend_pot = zeros(t_epoch * srate + length(t) - 1, 1); % Initialize dendritic potential
    
        % Generate post-synaptic potentials
        for k = 1:psn_num
            tau_s = 0.001 * (2 * rand + 2); % Synaptic timescale (seconds)
            f = (t.^0.13) .* exp(-(t + (lat_s * srate)) / (tau_s * srate));
            f = f / sum(f); % Normalize PSP
    
            % Generate spike train for the epoch
            spike_train = zeros(t_epoch * srate, 1);
    
            if is_high
                % High firing rate for the first part
                for step = 1:high_samples
                    spike_train(step) = rand <= (frate_high / srate);
                end
            end
    
            % Low firing rate for the remainder of the epoch
            start_low = high_samples + 1;
            for step = start_low:(t_epoch * srate)
                spike_train(step) = rand < (frate_low / srate);
            end
    
            % Random sign and magnitude for spike train
            spike_train_withMag = spike_train .* (rand(size(spike_train)) - 0.5) * 2; 
            dend_pot = dend_pot + conv(spike_train_withMag, f, 'full'); % Convolve spike train with PSP
        end
    
        % Store cumulative dendritic potential for this epoch
        start_idx = (epoch - 1) * interval * srate + 1;
        end_idx = epoch * interval * srate;
    
        if epoch == 1
            time_series(start_idx:end_idx) = cumsum(dend_pot(1:(end_idx - start_idx + 1)));
        else
            time_series(start_idx:end_idx) = cumsum(dend_pot(1:(end_idx - start_idx + 1))) + ...
                time_series((epoch - 1) * interval * srate);
        end
    end

    time_series_all(:,un) = ieeg_highpass(time_series, srate);

end

end
