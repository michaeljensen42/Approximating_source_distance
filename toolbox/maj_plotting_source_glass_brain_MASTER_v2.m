% load channels near dipole, and estimate point source of dipole
%clear
pt = 'DIS';
stim_chan = 'RZ1';

load('/Users/M215487/Documents/WORK/SEEG/SEEG_source/data/DIS/chopped/carla/DIS_ccep_098_099_div_carla_hp_post.mat', 'V', 't', 'srate') % load bipolar re-referenced, high passed data
load(['brains' filesep pt '_brain.mat'],'locs', 'lbls');
    
    lead_names = {'RN', 'RY'};
    lbls_short = lbls; lbls_short(isnan(locs(:,1))) = []; % channel names for CARLA
    locs_short = locs; locs_short(isnan(locs(:,1)),:) = []; % carla channel stuff
    RN_RY_idx_all = find(contains(lbls_short, lead_names));
    RN_RY_idx = find(contains(lbls_short, lead_names));
        RN_RY_idx(22:25) = []; % get rid of RY channels in different cortex as they exist (RY11-RY14)
        RN_RY_idx(10:11) = []; % get rid of RN channels in different cortex as they exist (RN11-RN12)
        chanNames = lbls_short(RN_RY_idx);
        chanLocs =  locs_short(RN_RY_idx,:);

    V_tmp = V(:,:,RN_RY_idx); % index into data for channels of interest
    V_RN_RY = squeeze(reshape(V_tmp, [], size(V_tmp,1), size(V_tmp,2)*size(V_tmp,3))); % reshape into 1 2D matrix (time x trials)

    tpts=find(and(t>0.015 ,t<=2)); % start paramterizing at 15 ms after stim until 2 seconds after stim
    t_win=t(tpts);

    V=V_RN_RY(tpts,:);

    [crp_parms, crp_projs]=CRP_method_fig6(V,t_win); % run CRP treating all channels as 1 large set

    mean_V = zeros(length(RN_RY_idx), srate/2+1);
    for ch = 1:length(RN_RY_idx)

        al_p(ch,:) = -1*crp_parms.al_p((ch-1)*12+1:ch*12);
        mean_al_p(ch) = -1*mean(crp_parms.al_p((ch-1)*12+1:ch*12));
        [~,p_value_noCorrected(ch)] = ttest(crp_parms.al_p((ch-1)*12+1:ch*12), 0,'tail','left'); % get significance of alpha primes using 1 sample t-test vs. 0, correct for number of comparisons
        p_value = p_value_noCorrected.*length(RN_RY_idx);
        mean_V(ch,:) = mean(V_RN_RY(0.5*srate:srate,(ch-1)*12+1:ch*12),2)';

    end

    sig_chanNames = chanNames(p_value<0.05);
    sig_mean_al_p = mean_al_p(p_value<0.05)';
    sig_chanLocs = chanLocs(p_value<0.05,:);

    clear crp* elec_descriptions bmat bvol CAR ch cortex* x y z V* stats t t_win theta k _opt phi

[source_position, k_opt, theta, phi] = estimate_voltage_source_v3(sig_chanLocs, sig_mean_al_p);

locs_tmp = sig_chanLocs;
lbls_tmp = sig_chanNames;

locs_tmp(end+1,:) = source_position;
lbls_tmp(end+1) = {'SRC'};

 
%% Plot in slices, red + markes the source, Panel C position of source

    slicethickness = 8;
    mr_clims=[0 1];
    side = 'R';

    load(['brains' filesep pt '_brain.mat']);

% Plot Labels

   [x_slice, y_slice, z_slice] = seegview_sliceplot(locs_tmp, bvol, x, y, z, slicethickness, mr_clims, side); 
    
    sv_label_add_numbers(locs_tmp, lbls_tmp, x_slice);
    sv_label_add_numbers(locs_tmp, lbls_tmp, y_slice);
    sv_label_add_numbers(locs_tmp, lbls_tmp, z_slice);


   
%% Panel C, time seris from each channel
% plot significant channels based on alpha prime values with their distance to the source above each trace alongside the alpha prime value of that channel

      carla_dists2source = vecnorm(source_position - sig_chanLocs,2,2);
                carlam = max(carla_dists2source);
                wts_carla = 2*carla_dists2source(:)/carlam;
                    wts_carla = (wts_carla - min(wts_carla))/max(wts_carla); %shift so min is zero
                    wc = 0;

      gap = 50; 
      figure

       for jj = 1:length(RN_RY_idx)

           if p_value(jj) < 0.05

            wc = wc+1;
            hold on, plot(gap*(jj-1)+mean_V(jj,:), 'Color', [0+wts_carla(wc) 0+wts_carla(wc) 0+wts_carla(wc)], 'LineWidth', 1);                            
            hold on, text(srate/5, gap*(jj-1), sprintf('%3.3f %3.3f', sig_mean_al_p(wc), carla_dists2source(wc)), 'FontSize', 10);
            clear p l1

           else 

            hold on, plot(gap*(jj-1)+ mean_V(jj,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 1);

           end
       end
                                title('high pass, notch carla')
                                clear jj

                                xticks(0:srate/10:0.5*srate); xticklabels({'0ms', '100', '200', '300', '400', '500'})
                                yticks(0:gap:length(chanNames)*gap); yticklabels(chanNames)

            
