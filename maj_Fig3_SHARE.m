function maj_Fig3_SHARE
% load data

clear
pt = 'BYG';
reref = 'dp';
stim_pair = [139 140];
rej_method = 'bathtub'; 

load(['brains' filesep pt '_brain' ],'locs','lbls');
[dp_channels, dp_locs]=locs_DPRR(locs); dp_lbls = lbls(dp_channels); 

chansOfI = {'RM2', 'RM4', 'RM6'};
ch_idx = find(contains(dp_lbls, chansOfI)');

load(['data' filesep pt filesep 'chopped' filesep reref filesep pt '_ccep_' sprintf('%3.3d',stim_pair(1)) '_' sprintf('%3.3d',stim_pair(2)) '_div_' reref '_epsilon_' rej_method])

load(['data' filesep pt filesep 'CRP' filesep reref filesep pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' sprintf('%3.3d',stim_pair(2)) '_div_' reref '_hp_post'])

load(['data' filesep pt filesep 'chopped' filesep reref filesep pt '_ccep_' sprintf('%3.3d',stim_pair(1)) '_' sprintf('%3.3d',stim_pair(2)) '_div_' reref '_hp_post'])

%% run wavlets and Plot panels B, E, H
figure

for ch = ch_idx  

    ch_name = dp_lbls{ch};
    frange = 1:200;
    method = 'kjm';
    clim_lower = 1;
    clim_upper = 1.7;

     [S_ep, f] = ieeg_getWaveletSpectrogram(squeeze(data_ep(:,ch)), srate, frange, method); % for 

        for st = 1:length(stim_indices)
            S_baseline(:,st,:) =  S_ep(:,stim_indices(st)-0.22*srate:stim_indices(st)-0.02*srate); % baseline -200 to 20 ms prior to stim
        end 
        clear st

        S_base2 = squeeze(mean(S_baseline,2));
        S_basevec = squeeze(mean(S_base2,2));
        S_norm = S_ep./S_basevec; % average entire run by 500 to 100 ms prior to stim baseline

        for sti = 1:length(stim_indices)
            S_chopped(:,sti,:) = S_norm(:,stim_indices(sti)-srate:stim_indices(sti)+srate);
        end
        clear sti


        time = (-20:1000/srate:500);
        S_avg_across_trials = squeeze(mean(S_chopped,2));

        spec_ts = S_avg_across_trials(:, 0.98*srate:1.5*srate); % XXX



            nexttile
            imagesc(time, [], spec_ts)
            hold on, 
            %xline(0.2, 'Color', 0.5*[1 1 1], 'LineWidth', 4)
            
            axis xy
            load dg_colormap.mat; cm = cm(33:end,:); colormap(cm);
            clim([clim_lower clim_upper]); colorbar

            ylabel('Frequency (Hz)'); xlabel('Time (s)'); set(gcf, 'Color', 'w');
            yticks(0:50:size(S_avg_across_trials,1));
            title(sprintf('%s %s', reref, ch_name))


end

%% plot parameterization, panels C, F, I, Bipolar
figure

for ch = ch_idx  

    nexttile

    V_nostim(:,bad_input_trials{ch}, :);

    meanV = mean(V_nostim(srate/2:1.5*srate, :, ch),2);

    l1 = plot(V_nostim(srate/2:srate, :, ch), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

    [~,p] = ttest(crp_parms(ch).al_p, 0,'tail','right');

    if p < 0.05
        hold on, l2 = plot(meanV(1:crp_parms(ch).tR*srate, 1), 'y', 'LineWidth', 6);
    end

    hold on, l3 = plot(meanV(1:srate/2), 'k', 'LineWidth', 3);

    hold on, l4 = plot([1 srate/2],[0 0],'color',.3*[1 1 1], 'LineWidth', 2);

        xticks([1, srate/4, srate/2]) % downsampled so srate/8 is 250 ms
        xticklabels({'15', '250', '500' });

        xlabel('time (ms)'); ylabel('µV');

            if strcmp(reref, 'dp')
    

                ylim([-250 250]);
                yticks([-100 0 100])
                yticklabels({'-100', '0', '100'})

                title(sprintf('%s to %s', lbls{stim_pair(1)}, dp_lbls{ch}))
    
            elseif strcmp(reref, 'carla') 
    
                % ylim([-250 250]);
                % yticks([-250 0 250])
                % yticklabels({'-250', '0', '250'})

                ylim([-250 250]);
                yticks([-100 0 100])
                yticklabels({'-100', '0', '100'})

                title(sprintf('%s to %s', lbls{stim_pair(1)}, lbls_short{ch}))
    
            end

        xlim([-15 srate/2+0.015*srate]);

            if p < 0.05
                legend([l1(1), l3 ,l2], 'all trials', 'average trace', 'tR of average trace'); legend('boxoff'); box off; set(gcf, 'Color', 'w')
            else
                legend([l1(1), l3], 'all trials', 'average trace', 'tR of average trace'); legend('boxoff'); box off; set(gcf, 'Color', 'w')
            end

        set(gcf, 'Color', 'w'); box off

end

%% load CARLA data

clear
pt = 'BYG';
reref = 'carla';
stim_pair = [139 140];
rej_method = 'bathtub'; 

load(['brains' filesep pt '_brain' ],'locs','lbls');
lbls_short = lbls; lbls_short(isnan(locs(:,1))) = []; % channel names for CARLA

chansOfI = {'RM2', 'RM4', 'RM6'};
ch_idx = find(contains(lbls_short, chansOfI)');

load(['data' filesep pt filesep 'chopped' filesep reref filesep pt '_ccep_' sprintf('%3.3d',stim_pair(1)) '_' sprintf('%3.3d',stim_pair(2)) '_div_' reref '_epsilon_' rej_method])

load(['data' filesep pt filesep 'CRP' filesep reref filesep pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' sprintf('%3.3d',stim_pair(2)) '_div_' reref '_hp_post'])

%% plot parameterization, panels D, G, J, CARLA

figure

for ch = ch_idx  

    nexttile

    V_nostim(:,bad_input_trials{ch}, :);

    meanV = mean(V_nostim(srate/2:1.5*srate, :, ch),2);

    l1 = plot(V_nostim(srate/2:srate, :, ch), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

    [~,p] = ttest(crp_parms(ch).al_p, 0,'tail','right');

    if p < 0.05
        hold on, l2 = plot(meanV(1:crp_parms(ch).tR*srate, 1), 'y', 'LineWidth', 6);
    end

    hold on, l3 = plot(meanV(1:srate/2), 'k', 'LineWidth', 3);

    hold on, l4 = plot([1 srate/2],[0 0],'color',.3*[1 1 1], 'LineWidth', 2);

        xticks([1, srate/4, srate/2]) % downsampled so srate/8 is 250 ms
        xticklabels({'15', '250', '500' });

        xlabel('time (ms)'); ylabel('µV');

            if strcmp(reref, 'dp')
    

                ylim([-250 250]);
                yticks([-100 0 100])
                yticklabels({'-100', '0', '100'})

                title(sprintf('%s to %s', lbls{stim_pair(1)}, dp_lbls{ch}))
    
            elseif strcmp(reref, 'carla') 
    
                ylim([-250 250]);
                yticks([-100 0 100])
                yticklabels({'-100', '0', '100'})

                title(sprintf('%s to %s', lbls{stim_pair(1)}, lbls_short{ch}))
    
            end

        xlim([-15 srate/2+0.015*srate]);

            if p < 0.05
                legend([l1(1), l3 ,l2], 'all trials', 'average trace', 'tR of average trace'); legend('boxoff'); box off; set(gcf, 'Color', 'w')
            else
                legend([l1(1), l3], 'all trials', 'average trace', 'tR of average trace'); legend('boxoff'); box off; set(gcf, 'Color', 'w')
            end

        set(gcf, 'Color', 'w'); box off

end

