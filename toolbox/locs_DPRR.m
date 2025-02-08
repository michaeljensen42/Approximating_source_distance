function [dp_channels, dp_locs]=locs_DPRR(locs)


%%
    sep_locs=sum(diff(locs).^2,2).^.5;
    dp_channels=find(sep_locs<5);

%%
    dp_locs=(locs(1:(end-1),:)+locs(2:end,:))/2;
    dp_locs=dp_locs(dp_channels,:);
    
%%
    %figure, plot(sep_locs,'ro'), hold on, plot(dp_channels,sep_locs(dp_channels,:),'go'),
    %xlabel('electrode difference pairs'), ylabel('distance between pair')
    %
    %figure, plot3(locs(:,1),locs(:,2),locs(:,3),'b.')
    %hold on, plot3(dp_locs(:,1),dp_locs(:,2),dp_locs(:,3),'r.')





