function [data_out] = maj_stim_artifact_rejection_v2(stim_indices, data, srate, method, buffer, reref)

    % inputs:
        % stim - double, samples x 1 array containg the timing of stimulation
        % data - double, samples x channels, high passed, notch filtered data prior to referencing
        % srate - sampling rate of experiment
        % buffer - buffer before and after stimulation that is to be diminshed using inverted hann window
            
    % output:
        % data_out = double, input data without stimulation artifact

    % method:
        % hann = apply an inverted hann window to the data centered at the stimulation pulse, the buffer determines the half-width of the hann
        % crowther = use the method outlines in Crowther et al

    if nargin < 5, buffer = 0.015; end % 15 msec before/after stim by default

    data_out = data; % make copy of data for clean up
        
    if strcmp(method, 'hann')

    % define the peri-stim indices to be cleaned up based on buffer size
    peri_stim_lims = [-0.8*buffer 1.2*buffer]; % define rejection window
    peri_stim_indices= floor(srate*peri_stim_lims(1)):floor(srate*peri_stim_lims(2)); 
                        
                        hann_inverted = 1-hann(2*buffer*srate+1);
            
                        for jj = 1:size(data,2) % loop through each channel
        
                            for k=1:length(stim_indices) % loop through each stimulation
                            
                                stim_data_tmp = data(stim_indices(k) + peri_stim_indices, jj); % data stim window for 1 chan, 1 stim
                
                                cleaned_data = hann_inverted.*stim_data_tmp; % multiply by inverted hann
        
                                data_out(stim_indices(k)+peri_stim_indices,jj) = cleaned_data; % stitch back in cleaned data
            
                                clear a b stim_data_tmp cleaned_data
        
                            end
        
                        end

    elseif strcmp(method, 'bathtub')

    % define the peri-stim indices to be cleaned up based on buffer size
    peri_stim_lims = [-0.9*buffer 2.1*buffer]; % define rejection window
    peri_stim_indices= floor(srate*peri_stim_lims(1)):floor(srate*peri_stim_lims(2)+1); 
                        
                        hann_inverted = 1-hann(2*buffer*srate+1);

                        %bathtub_tmp = [hann_inverted; hann_inverted];
                        %bathtub_tmp(buffer*srate+2:3*buffer*srate+2) = 0;

                        zero_vec = zeros(1.*buffer*srate,1);

                        bathtub = [hann_inverted(1:(buffer*srate)); zero_vec; hann_inverted(buffer*srate:length(hann_inverted))];
            
                        for jj = 1:size(data,2) % loop through each channel
        
                            for k=1:length(stim_indices) % loop through each stimulation
                            
                                stim_data_tmp = data(stim_indices(k) + peri_stim_indices, jj); % data stim window for 1 chan, 1 stim
                
                                cleaned_data = bathtub.*stim_data_tmp; % multiply by inverted hann
        
                                data_out(stim_indices(k)+peri_stim_indices,jj) = cleaned_data; % stitch back in cleaned data
            
                                clear a b stim_data_tmp cleaned_data
        
                            end
        
                        end


    elseif strcmp(method, 'crowther')


    pre_stim_lims = [-buffer 0]; % define rejection window
    pre_stim_indices=round(srate*pre_stim_lims(1)):round(srate*pre_stim_lims(2)); 
    pre_wt = 1:-1/(buffer*srate):0; % weight vector for pre stim data

    post_stim_lims = [buffer 2*buffer]; % define rejection window
    post_stim_indices=round(srate*post_stim_lims(1)):round(srate*post_stim_lims(2)); 
    post_wt = 0:1/(buffer*srate):1; % weight vector for post stim data

    peri_stim_lims = [0 1*buffer]; % define rejection window
    peri_stim_indices= floor(srate*peri_stim_lims(1)):floor(srate*peri_stim_lims(2)); 

                     for jj = 1:size(data,2) % loop through each channel
            
                            for k=1:length(stim_indices)
                                        
                                    pre_tmp = data(stim_indices(k) + pre_stim_indices, jj); % data from pre-stim window
                                    rpre_tmp = pre_tmp(end:-1:1); % flipping pre_stim window
                                    a = rpre_tmp.*pre_wt'; % weighting the flipped pre stim window 
                
                                    post_tmp = data(stim_indices(k) + post_stim_indices,jj); % as above for post stim window
                                    rpost_tmp = post_tmp(end:-1:1); 
                                    b = rpost_tmp.*post_wt';
                                    
                                    
                
                                    combo = a+b;
                                                
                                    data_out(stim_indices(k)+peri_stim_indices,jj) = combo';

                                                                        clear a b pre_tmp rpre* post_tmp rpost*
            
                            end
            
                     end

    elseif strcmp(method, 'pca')


    buffer = 0.004;

    peri_stim_lims = [-0.6*buffer 1.9*buffer]; % define rejection window
    peri_stim_indices= floor(srate*peri_stim_lims(1)):floor(srate*peri_stim_lims(2))-1; 

                     hann_win = hann(2*buffer*srate)';

                     ones_vec = ones(1,srate/100 - length(hann_win));

                     bathtub = [hann_win(1:(buffer*srate)), ones_vec, hann_win(buffer*srate:length(hann_win))];

                     hann_inverted = 1-hann(2*buffer*srate)';

                        zero_vec = zeros(1,srate/100 - length(hann_win));

                        inverted_bathtub = [hann_inverted(1:(buffer*srate)), zero_vec, hann_inverted(buffer*srate:length(hann_inverted))];

                     for jj = 1:size(data,2) % loop through each channel
        
                            for k=1:length(stim_indices) % loop through each stimulation
                            
                                if strcmp(reref, 'dp')
                                    stim_data_tmp(k,:) = data(stim_indices(k) + peri_stim_indices, jj) - mean(data(stim_indices(k)-0.02*srate:stim_indices(k)-0.01*srate, jj)); % data stim window for 1 chan, 1 stim
                                else
                                    stim_data_tmp(k,:) = data(stim_indices(k) + peri_stim_indices, jj); % data stim window for 1 chan, 1 stim
                                end
                
                                %hann_data(:,k) = hann_win.*stim_data_tmp(k,:); % multiply by inverted hann

                                hann_data(:,k) = bathtub.*stim_data_tmp(k,:); % multiply by inverted hann

        
                            end

                                [F,S2]=eig(hann_data.'*hann_data); % eigenvector decomposition of (covariance of transpose) - may want to replace this with cov function so mean is subtracted off
                                [S2,v_inds]=sort(sort(sum(S2)),'descend'); F=F(:,v_inds); %reshape properly    
                                %
                                S=S2.^.5; % estimated eigenvalues of both X.'*X and X*X.'

                                ES=hann_data*F; % kernel trick
                                Eigvecs =ES./(ones(size(hann_data,1),1)*S);
                                stim_pc1 = Eigvecs(:,1);
                            
                                wt=stim_pc1.'*hann_data; % alpha coefficient weights for C into V
                                ep= hann_data-stim_pc1*wt; % residual epsilon (error timeseries) after removal of form of CCEP 

                                data_mid = data; % create additional copy for middle step
                                        for k=1:length(stim_indices) % loop through each stimulation
            
                                            data_out(stim_indices(k) + peri_stim_indices, jj) = data(stim_indices(k) + peri_stim_indices, jj) - (wt(k)*stim_pc1); % stitch back in cleaned data

                                            %data_out(stim_indices(k) + peri_stim_indices, jj) = data_mid(stim_indices(k) + peri_stim_indices, jj).*inverted_bathtub'; % re-weight using inverted bath tub

                                        end

                                                                                    clear a b stim_data_tmp ep wt stim_pc1 ES S F S2 V_inds
                                    
                      end







    end

end
