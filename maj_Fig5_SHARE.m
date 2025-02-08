function maj_Fig5_SHARE

    srate = 4800; 

    % panel A

        pt = 'DIS';
        stim_pair = [186 187];
    
        figure, subplot(211)
    
            load(['data/' pt '/' 'CRP' '/carla/' pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' ...
                    sprintf('%3.3d',stim_pair(2)) '_div_carla_hp.mat']);
    
        hold on, plot(100 + mean(V(0.515*srate:1*srate,:, 25),2), 'k')
        hold on, plot(mean(V(0.515*srate:1*srate,:, 24),2), 'k')
        title('mCAR')
        ylim([-200 200])
        set(gcf, 'Color', 'w'); box off
    
        clear V
    
            load(['data/' pt '/' 'CRP' '/dp/' pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' ...
                sprintf('%3.3d',stim_pair(2)) '_div_dp_hp.mat']);
    
        subplot(212)
        hold on, plot(100 + mean(V(0.515*srate:1*srate,:, 21),2), 'k')
        title('bipolar')
        ylim([-200 200])
        set(gcf, 'Color', 'w'); box off

    % panel B

        clear
        srate = 4800; 
        pt = 'BYG';
        stim_pair = [138 139];
    
        figure, subplot(211)
    
            load(['data/' pt '/' 'CRP' '/carla/' pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' ...
                    sprintf('%3.3d',stim_pair(2)) '_div_carla_hp.mat']);
    
        hold on, plot(100 + mean(V(0.515*srate:1*srate,:, 45),2), 'k')
        hold on, plot(mean(V(0.515*srate:1*srate,:, 44),2), 'k')
        title('mCAR')
        ylim([-200 500])
        set(gcf, 'Color', 'w'); box off
    
        clear V
    
            load(['data/' pt '/' 'CRP' '/dp/' pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' ...
                sprintf('%3.3d',stim_pair(2)) '_div_dp_hp.mat']);
    
        subplot(212)
        hold on, plot(100 + mean(V(0.515*srate:1*srate,:, 34),2), 'k')
        title('bipolar')
        ylim([-200 500])
        set(gcf, 'Color', 'w'); box off

    % panel C

        clear
        srate = 4800; 
        pt = 'USB';
        stim_pair = [102 103];
    
        figure, subplot(211)
    
            load(['data/' pt '/' 'CRP' '/carla/' pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' ...
                    sprintf('%3.3d',stim_pair(2)) '_div_carla_hp.mat']);
    
        hold on, plot(100 + mean(V(0.515*srate:1*srate,:, 54),2), 'k')
        hold on, plot(mean(V(0.515*srate:1*srate,:, 53),2), 'k')
        title('mCAR')
        ylim([-200 500])
        set(gcf, 'Color', 'w'); box off
    
        clear V
    
            load(['data/' pt '/' 'CRP' '/dp/' pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' ...
                sprintf('%3.3d',stim_pair(2)) '_div_dp_hp.mat']);
    
        subplot(212)
        hold on, plot(100 + mean(V(0.515*srate:1*srate,:, 43),2), 'k')
        title('bipolar')
        ylim([-200 500])
        set(gcf, 'Color', 'w'); box off

    % panel D

        srate = 4800; 
        pt = 'BYG';
        stim_pair = [138 139];
    
        figure, subplot(211)
    
            load(['data/' pt '/' 'CRP' '/carla/' pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' ...
                    sprintf('%3.3d',stim_pair(2)) '_div_carla_hp.mat']);
    
        hold on, plot(100 + mean(V(0.515*srate:1*srate,:, 88),2), 'k')
        hold on, plot(mean(V(0.515*srate:1*srate,:, 87),2), 'k')
        title('mCAR')
        ylim([-200 500])
        set(gcf, 'Color', 'w'); box off
    
        clear V
    
            load(['data/' pt '/' 'CRP' '/dp/' pt '_CRP_' sprintf('%3.3d',stim_pair(1)) '_' ...
                sprintf('%3.3d',stim_pair(2)) '_div_dp_hp.mat']);
    
        subplot(212)
        hold on, plot(100 + mean(V(0.515*srate:1*srate,:, 72),2), 'k')
        title('bipolar')
        ylim([-200 500])
        set(gcf, 'Color', 'w'); box off

end