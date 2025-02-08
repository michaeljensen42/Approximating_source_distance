function maj_Fig4_SHARE()


% to generate new time series from scratch. NOTE: this will produce unique
% results from Figure 4
trials = 20;
units = 10;

constant_frate_low = 2; % for the noise


% generate noise using psd simulation code

[noise_time_series_low, srate] = generate_noise_time_series_v2(constant_frate_low, trials);
%[noise_time_series_med, ~] = generate_noise_time_series_v2(constant_frate_medium);
%[noise_time_series_high, ~] = generate_noise_time_series_v2(constant_frate_high);

multiplier = 3;
[whitenoise, ~] = generate_whitebrownnoise_time_series_v1(multiplier, trials); % instrument noise

% Broadband time series using psd simulation code and noise

% ERP time series and noise
amplitude_multiplier_high = 1.;
amplitude_multiplier_medium = 0.25;
amplitude_multiplier_low = 0.11;

frate_high = 2; % stim based firing
frate_base = 0; % stim based firing

[BB_time_series, srate, t_tot, interval] = generate_BB_time_series_v4(frate_high, frate_base, trials, units);
 
BB_time_series_high = BB_time_series * amplitude_multiplier_high;
BB_time_series_medium = BB_time_series * amplitude_multiplier_medium;
BB_time_series_low = BB_time_series * amplitude_multiplier_low;

erp_time_series = generate_erp_time_series_v2(trials, units); % generate ERP with same background firing rate, noise as BB_time_series

erp_time_series_high = erp_time_series * amplitude_multiplier_high;
erp_time_series_medium = erp_time_series * amplitude_multiplier_medium;
erp_time_series_low = erp_time_series * amplitude_multiplier_low;

compound_time_series_high = sum([BB_time_series_high, noise_time_series_low, erp_time_series_high, whitenoise], 2);
compound_time_series_med = sum([BB_time_series_medium, noise_time_series_low, erp_time_series_medium, whitenoise], 2);
compound_time_series_low = sum([BB_time_series_low, noise_time_series_low, erp_time_series_low, whitenoise], 2);


clear fr*

%% LOAD exact time series used to make figure 4 in manusrcipt
load('data/Fig4_timeseries.mat')


%% plot pure erp, noise, and BB times series for bottom left panel 4 TRIALs

stim_indices = [interval*srate:interval*srate:(t_tot-interval)*srate];

figure
plot(erp_time_series_high(stim_indices(1)-srate:stim_indices(4)+srate), 'k')
ylim([-200 100])
title('short erp ts')


figure
plot(BB_time_series_high(stim_indices(1)-srate:stim_indices(10)+srate), 'k')
ylim([-200 100])
title('short BB ts')

figure
plot(noise_time_series_low(stim_indices(1)-srate:stim_indices(4)+srate), 'k')
ylim([-200 100])
title('short noise ts')

figure
plot(whitenoise(stim_indices(1)-srate:10:stim_indices(4)+srate), 'k')
ylim([-200 100])
title('short instrument noise ts')



%% plot compound signals with low high medium SNR FULL TIME

figure
plot(compound_time_series_low(1:100:end), 'k')
title('full low SNR')
ylim([-300 200])


figure
plot(compound_time_series_med(1:100:end), 'k')
title('full med SNR')
ylim([-300 200])



figure
plot(compound_time_series_high(1:100:end), 'k')
title('full high SNR')
ylim([-300 200])


        %% Run CRP on erp_time_series and plot, high SNR

        stim_indices = interval*srate:interval*srate:(t_tot-interval)*srate;

        compound_time_series_high_spec = maj_stim_artifact_rejection_v2(stim_indices, compound_time_series_high, srate, 'bathtub');

                    for kk = 1:length(stim_indices)
                        % stack trial with baseline correction
                        V_plot(:,kk) = compound_time_series_high(stim_indices(kk):stim_indices(kk)+srate-1); % - mean(compound_time_series_high(stim_indices(kk)-0.05*srate:stim_indices(kk))); %mean(erp_time_series_low(stim_indices(kk)-srate/5:stim_indices(kk)-srate/10));
                        V_spec(:,kk) = compound_time_series_high_spec(stim_indices(kk):stim_indices(kk)+srate-1); % - mean(compound_time_series_high_spec(stim_indices(kk)-0.05*srate:stim_indices(kk))); %mean(erp_time_series_low(stim_indices(kk)-srate/5:stim_indices(kk)-srate/10));
                    end
                   
                    t_win = (1:size(V_spec,1))/srate;
            
                   [crp_parms_plot, crp_projs_plot] = CRP_method(V_plot,t_win);
                   [crp_parms_spec, crp_projs_spec] = CRP_method(V_spec,t_win);
            
                   figure

                   meanV = mean(V_plot,2);

                    l1 = plot(V_plot(1:10:srate/2,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

                    hold on, l2 = plot(crp_projs_plot.avg_trace_input(1:10:crp_parms_plot.tR*srate), 'y', 'LineWidth', 8);

                    hold on, l3 = plot(crp_projs_plot.avg_trace_input(1:10:srate/2), 'k', 'LineWidth', 3);

                    hold on, l4 = plot([1 srate/20],[0 0],'color',.3*[1 1 1], 'LineWidth', 2);

                    xticks([1, srate/4, srate/2]) % downsampled so srate/8 is 250 ms
                    xticklabels({'0', '250', '500' });

                    xlabel('time (ms)'); ylabel('µV');

                            ylim([-350 250]);
                            yticks([-250 0 250])
                            yticklabels({'-250', '0', '250'})

                            title('Noiseless ERP time series low noise')

                    %xlim([-15 srate/2+0.015*srate]);
                    legend([l1(1), l3 ,l2], 'all trials', 'average trace', 'tR of average trace'); 
                    legend('boxoff'); box off; set(gcf, 'Color', 'w')


 
           stim_indices = interval*srate:interval*srate:(t_tot-interval)*srate;
            
           data_highSNR_with_epsilon1 = compound_time_series_high_spec; % make copy stitch into (remove synchronous activity)
           data_highSNR_with_epsilon2 = compound_time_series_high; % make copy stitch into (remove synchronous activity)

           % stitch in epsilon
           for jj = 1:length(stim_indices)
                
                % start_data = stim_indices(jj) + round(0.015*srate);
                % stop_data = stim_indices(jj) + round(0.015*srate) + length(crp_parms.ep(:,jj))-1;

                start_data = stim_indices(jj);
                stop_data = stim_indices(jj)+ size(crp_parms_spec.ep,1)-1;
    
                data_highSNR_with_epsilon1(start_data:stop_data) = crp_parms_spec.ep(:,jj);
                data_highSNR_with_epsilon2(start_data:stop_data) = data_highSNR_with_epsilon2(start_data:stop_data) - (crp_parms_spec.al(jj)*crp_parms_spec.C); 
    
                clear start* stop*
        
            end
       
     %  Generate and plot spectrograms, high SNR

        frange = 1:200;
        method = 'kjm';

         [S_ep, f] = ieeg_getWaveletSpectrogram(data_highSNR_with_epsilon1, srate, frange, method);
        %[Sn, f] = ieeg_getWaveletSpectrogram(data_highSNR_with_epsilon2, srate, frange, method);
        %[Sn, f] = ieeg_getWaveletSpectrogram(compound_time_series_high, srate, frange, method);

            for st = 1:length(stim_indices)
                S_baseline(:,st,:) =  S_ep(:,stim_indices(st)-0.22*srate:stim_indices(st)-0.02*srate); % baseline -200 to 20 ms prior to stim
            end 


            S_base2 = squeeze(mean(S_baseline,2));
            S_basevec = squeeze(mean(S_base2,2));
            S_norm = S_ep./S_basevec; % average entire run by 500 to 100 ms prior to stim baseline
        
        %Sn_mean = mean(Sn, 2);

        %Sn_norm = Sn./Sn_mean;

        stim_indices = [interval*srate:interval*srate:(t_tot-interval)*srate];
        for kk = 1:length(stim_indices)
            Sn_chopped(kk,:,:) = S_norm(:,stim_indices(kk)-srate/5:stim_indices(kk)+srate/5);
        end

            Sn_avg_trial = squeeze(mean(Sn_chopped));

            figure

            imagesc(Sn_avg_trial)
            axis xy
            load dg_colormap.mat; cm = cm(33:end,:); colormap(cm);
            clim([1 1.7]); colorbar

            xticks(1:0.2*srate:0.4*srate+1)
            xticklabels({'-200', '0', '200'})

            ylabel('Frequency (Hz)'); xlabel('Time (ms)'); set(gcf, 'Color', 'w');


        clear S* stim_i* method kk jj
        clear l* mean* t_win crp* kk V data*





        %% Run CRP on erp_time_series and plot, medium SNR

        

                stim_indices = [interval*srate:interval*srate:(t_tot-interval)*srate];

                compound_time_series_med_spec = maj_stim_artifact_rejection_v2(stim_indices, compound_time_series_med, srate, 'bathtub');

                for kk = 1:length(stim_indices)
                    % stack trial with baseline correction
                    V_plot(:,kk) = compound_time_series_med(stim_indices(kk):stim_indices(kk)+srate-1); 
                    V_spec(:,kk) = compound_time_series_med_spec(stim_indices(kk):stim_indices(kk)+srate-1); 
                end
               
                t_win = (1:size(V_spec,1))/srate;
        
                   [crp_parms_plot, crp_projs_plot] = CRP_method(V_plot,t_win);
                   [crp_parms_spec, crp_projs_spec] = CRP_method(V_spec,t_win);
        
        
               figure
        
                meanV = mean(V_plot,2);
        
            l1 = plot(V_plot(1:10:srate/2,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        
            hold on, l2 = plot(crp_projs_plot.avg_trace_input(1:10:crp_parms_plot.tR*srate), 'y', 'LineWidth', 6);
        
            hold on, l3 = plot(crp_projs_plot.avg_trace_input(1:10:srate/2), 'k', 'LineWidth', 3);
        
            hold on, l4 = plot([1 srate/20],[0 0],'color',.3*[1 1 1], 'LineWidth', 2);
        
                xticks([1, srate/4, srate/2]) % downsampled so srate/8 is 250 ms
                xticklabels({'0', '250', '500' });
        
                xlabel('time (ms)'); ylabel('µV');
            
                        ylim([-250 250]);
                        yticks([-250 0 250])
                        yticklabels({'-250', '0', '250'})
        
                        title('Noiseless ERP time series low noise')
        
                %xlim([-15 srate/2+0.015*srate]);
                legend([l1(1), l3 ,l2], 'all trials', 'average trace', 'tR of average trace'); 
                legend('boxoff'); box off; set(gcf, 'Color', 'w')
        

             
           data_medSNR_with_epsilon1 = compound_time_series_med_spec; % make copy stitch into (remove synchronous activity)
           data_medSNR_with_epsilon2 = compound_time_series_med; % make copy stitch into (remove synchronous activity)

           % stitch in epsilon
           for jj = 1:length(stim_indices)
                
                start_data = stim_indices(jj);
                stop_data = stim_indices(jj)+ size(crp_parms_spec.ep,1)-1;
    
                data_medSNR_with_epsilon1(start_data:stop_data) = crp_parms_spec.ep(:,jj);
                %data_medSNR_with_epsilon2(start_data:stop_data) = data_medSNR_with_epsilon2(start_data:stop_data) - (crp_parms.al(jj)*crp_parms.C); 
    
                clear start* stop*
        
            end
       
        % spectrograms, high SNR

        frange = 1:200;
        method = 'kjm';

        [S_ep, f] = ieeg_getWaveletSpectrogram(data_medSNR_with_epsilon1, srate, frange, method);
        %[Sn, f] = ieeg_getWaveletSpectrogram(data_medSNR_with_epsilon2, srate, frange, method);
        %[Sn, f] = ieeg_getWaveletSpectrogram(compound_time_series_med, srate, frange, method);

          for st = 1:length(stim_indices)
                S_baseline(:,st,:) =  S_ep(:,stim_indices(st)-0.22*srate:stim_indices(st)-0.02*srate); % baseline -200 to 20 ms prior to stim
          end 

            S_base2 = squeeze(mean(S_baseline,2));
            S_basevec = squeeze(mean(S_base2,2));
            S_norm = S_ep./S_basevec; % average entire run by 500 to 100 ms prior to stim baseline
        
        
        %Sn_mean = mean(Sn, 2);
        %Sn_norm = Sn./Sn_mean;

        stim_indices = [interval*srate:interval*srate:(t_tot-interval)*srate];
        for kk = 1:length(stim_indices)
            Sn_chopped(kk,:,:) = S_norm(:,stim_indices(kk)-srate/5:stim_indices(kk)+srate/5);
        end

            Sn_avg_trial = squeeze(mean(Sn_chopped));

            figure
            imagesc(Sn_avg_trial)
            axis xy
            
            %load dg_colormap.mat; cm = cm(33:end,:); 
            colormap(cm);
            clim([1 1.7]); colorbar

            xticks(1:0.2*srate:0.4*srate+1)
            xticklabels({'-200', '0', '200'})

            ylabel('Frequency (Hz)'); xlabel('Time (ms)'); set(gcf, 'Color', 'w');

        clear S* stim_i* method kk 

        clear l* mean* t_win crp* kk V

        %% Run CRP on erp_time_series and plot, low SNR

        stim_indices = [interval*srate:interval*srate:(t_tot-interval)*srate];

        compound_time_series_low_spec = maj_stim_artifact_rejection_v2(stim_indices, compound_time_series_low, srate, 'bathtub');

        for kk = 1:length(stim_indices)
            % stack trial with baseline correction
             V_plot(:,kk) = compound_time_series_low(stim_indices(kk):stim_indices(kk)+srate-1); 
             V_spec(:,kk) = compound_time_series_low_spec(stim_indices(kk):stim_indices(kk)+srate-1); 
        end
       
                t_win = (1:size(V_spec,1))/srate;
        
                   [crp_parms_plot, crp_projs_plot] = CRP_method(V_plot,t_win);
                   [crp_parms_spec, crp_projs_spec] = CRP_method(V_spec,t_win);


       figure

        meanV = mean(V_plot,2);

    l1 = plot(V_plot(1:10:srate/2,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

    hold on, l2 = plot(crp_projs_plot.avg_trace_input(1:10:crp_parms_plot.tR*srate), 'y', 'LineWidth', 6);

    hold on, l3 = plot(crp_projs_plot.avg_trace_input(1:10:srate/2), 'k', 'LineWidth', 3);

    hold on, l4 = plot([1 srate/20],[0 0],'color',.3*[1 1 1], 'LineWidth', 2);

        xticks([1, srate/4, srate/2]) % downsampled so srate/8 is 250 ms
        xticklabels({'0', '250', '500' });

        xlabel('time (ms)'); ylabel('µV');
    
                ylim([-250 250]);
                yticks([-250 0 250])
                yticklabels({'-250', '0', '250'})

                title('Noiseless ERP time series low noise')

        %xlim([-15 srate/2+0.015*srate]);
        legend([l1(1), l3 ,l2], 'all trials', 'average trace', 'tR of average trace'); 
        legend('boxoff'); box off; set(gcf, 'Color', 'w')


           data_lowSNR_with_epsilon1 = compound_time_series_low_spec; % make copy stitch into (remove synchronous activity)
           data_lowSNR_with_epsilon2 = compound_time_series_low; % make copy stitch into (remove synchronous activity)

           % stitch in epsilon
           for jj = 1:length(stim_indices)
                
                start_data = stim_indices(jj);
                stop_data = stim_indices(jj)+ size(crp_parms_spec.ep,1)-1;
    
                data_lowSNR_with_epsilon1(start_data:stop_data) = crp_parms_spec.ep(:,jj);
%                data_lowSNR_with_epsilon2(start_data:stop_data) = data_lowSNR_with_epsilon2(start_data:stop_data) - (crp_parms.al(jj)*crp_parms.C); 
    
                clear start* stop*
        
            end
       
        % spectrograms, high SNR

        frange = 1:200;
        method = 'kjm';

        [S_ep, f] = ieeg_getWaveletSpectrogram(data_lowSNR_with_epsilon1, srate, frange, method);
        %[Sn, f] = ieeg_getWaveletSpectrogram(data_lowSNR_with_epsilon2, srate, frange, method);
        %[Sn, f] = ieeg_getWaveletSpectrogram(compound_time_series_low, srate, frange, method);  

          for st = 1:length(stim_indices)
                S_baseline(:,st,:) =  S_ep(:,stim_indices(st)-0.22*srate:stim_indices(st)-0.02*srate); % baseline -200 to 20 ms prior to stim
            end 


            S_base2 = squeeze(mean(S_baseline,2));
            S_basevec = squeeze(mean(S_base2,2));
            S_norm = S_ep./S_basevec; % average entire run by 500 to 100 ms prior to stim baseline
        
        
        %Sn_mean = mean(Sn, 2);

        %Sn_norm = Sn./Sn_mean;

        stim_indices = [interval*srate:interval*srate:(t_tot-interval)*srate];
        for kk = 1:length(stim_indices)
            Sn_chopped(kk,:,:) = S_norm(:,stim_indices(kk)-0.2*srate:stim_indices(kk)+0.2*srate);
        end

            Sn_avg_trial = squeeze(mean(Sn_chopped));

            figure
            imagesc(Sn_avg_trial)
            axis xy
            %load dg_colormap.mat; cm = cm(33:end,:); 
            colormap(cm);
            clim([1 1.7]); colorbar

            xticks(1:0.2*srate:0.4*srate+1)
            xticklabels({'-200', '0', '200'})

            title('ep');

            ylabel('Frequency (Hz)'); xlabel('Time (ms)'); set(gcf, 'Color', 'w');



        clear S* stim_i* method kk 

        clear l* mean* t_win crp* kk V


end