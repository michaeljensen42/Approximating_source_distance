function sv_label_add_numbers(locs, loc_names, slice, msize)

    if nargin < 4, msize = 12; end % default marker (and text size)
    
    for k = 1:length(slice.values)
        axes(slice.ha(k))
        tmp_inds = slice.locs==k; % indices of electrodes corresponding to current slice
        sv_label_ax(locs(tmp_inds,slice.perm), loc_names(tmp_inds), msize)
    end
    
end

function sv_label_ax(els, lbls, msize)

    hold on, 

    scatterK.MarkerFaceAlpha = 0.9;
    for k=1:size(els,1) %k=1:2:length(els(:,1)) % text size 1/5th of marker size
        text(els(k,1), els(k,2), els(k,3), lbls{k}, 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'w')
            if strcmp(lbls{k}, 'SRC')
                scatterK = scatter3(els(k,1), els(k,2), els(k,3), msize, [1 0 0], 'LineWidth', 20);
            else
                scatterK = scatter3(els(k,1), els(k,2), els(k,3), msize, [0 1 0]);
            end
    end
end


        
        
        