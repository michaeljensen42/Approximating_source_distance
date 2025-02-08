function maj_Fig1_SHARE

load('data/USB/CRP/dp/USB_CRP_160_161_div_dp_hp_post.mat')

recording_chan = 70;

figure,

nexttile
    plot(crp_parms(recording_chan).V_tR)
nexttile
    plot(crp_parms(recording_chan).V_tR(:,3),'k')
    hold on, plot(crp_parms(recording_chan).ep(:,3), 'g', 'LineWidth', 2)
    hold on, plot(crp_parms(recording_chan).al(3)*crp_parms(recording_chan).C, 'r')
nexttile
    histogram(crp_parms(recording_chan).al_p, 'FaceColor', [0.5 0.5 0.5], 'BinWidth', 10)
    xlim([-10 130])

end