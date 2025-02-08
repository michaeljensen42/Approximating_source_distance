function [time_series, srate] = generate_noise_time_series_v2(frate, trials)

% input parameters
            %parameters
            srate=10000; %sampling rate (10kHz best?)
            lat_s=.002; %latency before decay kicks in (units seconds)
            t_tot= (trials*3) + 3; %total time (units of seconds)
            psn_num=6000; %number of presynaptic inputs (roughly 6,000-10,000 for a pyramidal neuron)
            t=0:floor(.05*srate); %50 ms window for PSP
            tau_s=.001*(2*rand+2);
            f=(t.^.13).*exp(-(t+(lat_s*srate))/(tau_s*srate));
            f=f/sum(f); 
            dend_pot=zeros(t_tot*srate+length(f)-1,1); %initialize dendritic potential
            
            for k=1:psn_num
                tau_s=.001*(2*rand+2); %synaptic timescale (units seconds)
                f=(t.^.13).*exp(-(t+(lat_s*srate))/(tau_s*srate)); f=f/sum(f); %create post-synaptic potential/current trace
                spike_times=(rand(t_tot*srate,1)<=(frate/srate)); %spike times - think of "a" as the trace at a single synapse
                spike_times_mag=spike_times.*(rand(t_tot*srate,1)-.5)*2; %arbitrary sign and magnitude on -1 to 1 for subsequent PSP/PSC
                dend_pot=dend_pot+conv(spike_times_mag,f); %convolve spike arrival  times with PSP shape - superimpose all synapses (i.e. passive dendrite approx.)
            end

             time_series_long =cumsum(dend_pot); %temporal integration - cumulative sum
             time_series = time_series_long(1:length(spike_times_mag))';

             time_series = ieeg_highpass(time_series, srate);

% Plot the generated time series
% t_plot = (0:length(time_series) - 1) / srate; % Time vector in seconds
% figure;
% plot(t_plot, time_series, 'k');
% xline(3:interval:t_tot);
% xlabel('Time (s)');
% ylabel('Cumulative Dendritic Potential');
% title('Asynchronous Time Series with Transitions Every 3 Seconds');


end