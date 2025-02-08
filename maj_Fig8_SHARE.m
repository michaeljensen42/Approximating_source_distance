function maj_Fig8_SHARE

    pt = 'RSM';
    rerefs = {'dp', 'carla', 'no_reref'};
    srate = 4800;
    stim_pair = [197 198];

    figure
        
    for ref = 1:length(rerefs)
       
            reref = rerefs{ref};

            load('data/dg_colormap.mat');
            
            if strcmp(reref, 'dp')

                load(['data/RSM/CRP/' reref '/RSM_CRP_197_198_div_' reref '_hp_post.mat'])
                load(['brains' filesep  pt '_brain'],'lbls', 'locs');[dp_channels,dp_locs]=locs_DPRR(locs); dp_lbls = lbls(dp_channels); % bipolar chan
                V(:,:, [60:61, 151, 155, 171, 10, 12, 20:21, 180, 181]) = []; % get rid of bad channels
            
            elseif strcmp(reref, 'carla')

                load(['data/RSM/CRP/' reref '/RSM_CRP_197_198_div_' reref '_hp_post.mat'])
                load(['brains' filesep  pt '_brain'], 'lbls', 'locs'); lbls_short = lbls; lbls_short(isnan(locs(:,1))) = [];
                 V(:,:, [73, 176, 181, 184]) = []; % get rid of bad channels

            elseif strcmp(reref, 'no_reref')

                load(['data/RSM/chopped/' reref '/RSM_ccep_197_198_div_' reref '_hp_post.mat'])
                load(['brains' filesep  pt '_brain'], 'lbls', 'locs');
                 V(:,:, [77, 84]) = []; % get rid of bad channels
            
            end
            
            a = squeeze(mean(V,2));
                
                nexttile
                imagesc(a(srate/4:srate,:)')
                colormap(cm)
                clim([-100 100])
            
                xticks(1:srate/4:size(V,1))
                xticklabels({'-250', '0', '250' ,'500'})
            
                xlabel('time (ms)', 'FontSize', 30);
                ylabel('channels', 'FontSize', 30)

                 if strcmp(reref, 'dp')
                    
                     title('bipolar')

                elseif strcmp(reref, 'carla')

                    title('carla')
              
                 elseif strcmp(reref, 'no_reref')

                    title('no re-reference')

                 end


    end     


end