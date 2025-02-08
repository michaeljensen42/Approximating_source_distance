function [erp_time_series] = generate_erp_time_series_v2(trials, units)

    srate = 10000;
    t_tot= (trials*3)+3;

    interval = 3;
    epochs = t_tot / interval; % Total number of epochs

    ERP_length = 0.2; % 200ms
    ERP_time = ERP_length*srate;
   

    time_series = zeros(t_tot*srate, 1); % Preallocate timeseries

 
        % Generate ERP
        dt = .001;
        t = 0:dt:2;
        x = linspace(1, length(t), length(t));
        
        k = 0.5;                   % Decay constant (adjustable)

        % Exponential curve equation
        %exp_decay = exp(-k * t) - exp(-2 * k);

        %gaus2 = [window(@gausswin,1.*floor(length(t)),3)]; %generate gaussian USE THIS
        gaus1 = [window(@gausswin,1.3*floor(length(t)),3)]; %generate gaussian USE THIS
        %gaus = [window(@gausswin,3*floor(length(t)),0.1)]; %generate gaussian
        %figure, plot(gaus(0.3*length(gaus):length(gaus)));
        %figure, plot(gaus);
        %gaushalf = gaus(end-length(gaus)/2+1:end); %cut gaussian in a half
        gauswin1 = gaus1(end-length(gaus1)/1.3:end); %first third to end
        %gauswin2 = gaus2(1:end); %first third to end

for un = 1:units

    % insert ERP into time series
    %figure
    for  epoch = 1:epochs-1
        % for each trial change amplitude and phase a bit
        % coeff1 = 1; %0.98 + (0.04* abs(randn));
        % coeff2 = 1; %0.98 + (0.04* abs(randn));

        sinwv1 = 15*sin(pi*1.5*t);
        sinwv2 = 3*sin(3*pi*1-t); %USE THIS
        %sinwv2 = 5*sin(3*pi*1-0.1*t);
    
        %ERP = exp_decay.*sinwv1.*sinwv2;
        ERP = gauswin1'.* sinwv1.* sinwv2; % .*gauswin2'
        %hold on, plot(ERP)
            
        start_time = epoch*srate*interval;
        time_series(start_time:start_time+ERP_time) = ERP;

        clear sin* coeff*
        
    end


                erp_time_series(:, un) = time_series;

end



% Plot the noise time series
% t_plot = (0:length(time_series) - 1) / srate; % Time vector in seconds
% figure;
% plot(t_plot, time_series, 'k');
% xlabel('Time (s)');
% ylabel('ERP Time Series');
% title('ERP Time Series without background firing');


end