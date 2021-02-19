% Pipeline compilation from 01/02/2021
% by Anton Tokariev (Baby Brain Activity Center/University of Helsinki)
%
% Original paper:
%'In utero exposure to antiepileptic drugs impacts sleep dynamics and
% developmental outcomes in infants'
%
% This pipeline computes frequency-specific phase-phase connectivity (PPC)
% for cortical signals from 19-channel EEG data
%
% Input: define the input EEG file (which is in .mat format)
%        see 'INPUT EEG DATA' section below
%        Example: load('Data example\filename.mat');]
%
% Output: PPC connectivity matrices computed using weighted phase lag index
% (for details on wPLI see: Vinck et al., 2011, Neuroimage)
%
% Copyright (C) 2021
% 
% This is a free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation
% 
% This program is distributed WITHOUT ANY WARRANTY
% (without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE)
%
% See the GNU General Public License: http://www.gnu.org/licenses


function compute_PPC()

 tic
      
 close all 
 clearvars

% INPUT EEG DATA
% ======================================================================= %
% eeg_data.eeg    - input EEG [channels x samples]
% eeg_data.Fs     - sampling rate
% eeg_data.Labels - channel names
  load('Data example\infant_eeg.mat'); %#ok<LOAD> 
 
  
% Visualization (of input EEG data):  
  plot_signals(eeg_data.eeg', eeg_data.Labels, eeg_data.Fs, 'EEG example', 'eeg');               
      
     
% LOAD FILTERS
% NOTE! Filters were computed for sampling rate Fs = 100 Hz!
% ======================================================================= %
% EEG filters (top: high-pass filters; bottom: low-pass filters)
% band 1: low delta  -> [0.4 - 1.5] Hz
% band 2: high delta -> [1.5 - 4] Hz
% band 3: theta      -> [4 - 8] Hz       
% band 4: alpha      -> [8 - 13] Hz
% band 5: beta       -> [13 - 22] Hz
  load('Filters\filters_eeg_5_bands.mat');   %#ok<LOAD>
    
  N_flt_eeg = size(flt_eeg, 2); % N of filters
  
  fr_bands = [{'low delta'}, {'high delta'}, {'theta'}, {'alpha'}, {'beta'}];

  
% LOAD HEAD MODEL
% ======================================================================= %         
  load('Head Model\Atlas.mat');                                            %#ok<LOAD>
  load('Head Model\InverseOperator.mat');                                  %#ok<LOAD>
  load('Head Model\CollapseOperator.mat');                                 %#ok<LOAD>
   
 
% ANALYSIS
% ======================================================================= %
% 1. Filter EEG
% ---------------------------------------------------
  disp('Filter EEG...');
  eeg_filtered = filter_eeg(eeg_data.eeg', flt_eeg);
% Output > eeg_filtered: cell{1 x N filters}(samples x channels)
  
% Visualization: filtered EEG
   fr = 3; % select number of frequency band
   plot_signals(eeg_filtered{1, fr}, eeg_data.Labels, eeg_data.Fs, ['Filtered EEG: ' fr_bands{1, fr}], 'eeg');
    
   
% 2. Compute parcel signals (see Tokariev et al., 2019, Cerebral Cortex)
% ---------------------------------------------------
% Visualization: show head model
  visualize_head_model()
  
  disp('Compute cortical signals...');
  parcels = get_parcel_signals(eeg_filtered, InverseOperator, CollapseOperator, Atlas);
% Output > parcels: cell{1 x N filters}(samples x parcels)

% Visualization: parcel signals
  fr = 3; % select number of frequency band
  plot_signals(parcels{1, fr}(:, 1:29),  Atlas.Areas(1:29),  eeg_data.Fs, ['Parcel signals/Left hemi: ' fr_bands{1, fr}],  'src');
  plot_signals(parcels{1, fr}(:, 30:58), Atlas.Areas(30:58), eeg_data.Fs, ['Parcel signals/Right hemi: ' fr_bands{1, fr}], 'src'); 
      
% 3. Phase-phase correlation (PPC) (see Vinck et al., 2011, Neuroimage)  
% ---------------------------------------------------
  disp('Compute PPC matrix...');
    
  PPC{1, N_flt_eeg} = [];                                                                              
  for n = 1:N_flt_eeg 
      PPC{1, n} = get_wPLI(parcels{1, n}); % PPC for each fr.band
      PPC{1, n} = PPC{1, n} + PPC{1, n}';  % make symmetric
  end
% Output > PPC: cell{1 x N filters}(parcels x parcels)
 
  PPC = cellfun(@abs, PPC, 'UniformOutput', false); % take abs 
   
% Visualization: PPC connectivity matrix
  fr = 3;
  figure; imagesc(PPC{1, fr}); colorbar; hold on; set(gcf, 'Color', 'w'); xlabel('parcels'); ylabel('parcels'); title(['PPC matrix: ' fr_bands{1, fr}], 'Fontweight', 'bold');  
   
  
% 4. PPC matrix correction (see Tokariev et al., 2019, Cerebral Cortex)
% ---------------------------------------------------
  disp('Correction of PPC matrix...');

% FidelityOperator: mask that removes 'non-reliable' edges  
  load('Head Model\FidelityOperator.mat');                                 %#ok<LOAD>
  
% Correct PPC matrix 
  for n = 1:N_flt_eeg
      PPC{1, n} = PPC{1, n} .* FidelityOperator;
  end
  
% Visualization: Corrected PPC connectivity matrix
  fr = 3;
  figure; imagesc(PPC{1, fr}); colorbar; c_map = get(gca, 'Colormap'); c_map(1, 1:3) = 0; set(gca, 'Colormap', c_map); hold on; set(gcf, 'Color', 'w'); xlabel('parcels'); ylabel('parcels'); title(['Corrected PPC matrix: ' fr_bands{1, fr}], 'Fontweight', 'bold');  
   
  
  
% EXAMPLE of 3D network Visualization
% ___________________________________
% Select frequency band
  fr = 3;
% Take strongest 5% of network connections (binary matrix)   
  PPC_top5 = double(PPC{1, fr} > prctile(nonzeros(PPC{1, fr}), 95));
% Visualization as (binary) connectivity matrix  
  figure; imagesc(PPC_top5); colormap(gray); colorbar; hold on; set(gcf, 'Color', 'w'); xlabel('parcels'); ylabel('parcels'); title(['Strongest 5% connections in ' fr_bands{1, fr} ' network'], 'Fontweight', 'bold');  
% Visualization as 3D network 
  show_3D_network(PPC_top5, Atlas, [fr_bands{1, fr} ' network (strongest 5%)']);
% ___________________________________  
  
% 5. Save PPC matrix  
% ---------------------------------------------------  
  disp('Save PPC matrix to "PPC output" folder...');
  
  save('PPC_output\PPC.mat', 'PPC');

 toc 

end % end of main function


function plot_signals(s, lables, Fs, ttl, type)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
%
% Input:
%      s: signal        [channels x samples]
% labels: channel names [channels x 1]
%     Fs: sampling rate
%    ttl: title (char, example ttl = 'Title') 
%   type: visualization  of EEG ('eeg') or parcel signals ('src')   

    L = size(s, 1); % length
   ch = size(s, 2); % number of channels
        
   shift = 5*(max(std(s))); % b/w different channels for visualization
   
   baseline = (1:1:ch)*shift;
   baseline = fliplr(repmat(baseline, L, 1));
  
 % time   
   t =(0:1:L-1)/Fs;
   t = repmat(t, ch, 1);
   t = t';
     
% add baseline to each channel for visualization
  sig_plot = s + baseline;
  
% plot signals
  figure;
  set(gcf, 'Color', 'w');
  hold on
  
% Plot EEG/Filtered EEG (all black) 
  if strcmp(type, 'eeg') == 1
     for i = 1:ch
         plot(t(:,i), sig_plot(:,i), 'Color', 'k');
     end
  end
  
% Plot parcel signals (colorful)   
  if strcmp(type, 'src') == 1
   % assign parcel colors   
     clr = zeros(ch, 3); % Occipital = black (default)
     clr(cell2mat(lables) == 'T', 1:3) = repmat([0.47 0.67 0.19],sum(cell2mat(lables) == 'T'), 1); % Temporal
     clr(cell2mat(lables) == 'C', 1:3) = repmat([0.75 0.00 0.75],sum(cell2mat(lables) == 'C'), 1); % Central
     clr(cell2mat(lables) == 'F', 1:3) = repmat([0.93 0.69 0.13],sum(cell2mat(lables) == 'F'), 1); % Frontal
     
     for i = 1:ch
         plot(t(:,i), sig_plot(:,i), 'Color', clr(i, :));
     end
  end
  
  xlabel('Time, sec', 'Fontweight', 'bold');
  ylabel('Channels',  'Fontweight', 'bold');
  
  set(gca, 'YTick', fliplr(baseline(1,:)), 'YTickLabel', flipud(lables));
  
  title(ttl); % from input
  
  Y_max = max(max(sig_plot));
  Y_min = min(min(sig_plot));  
  
  set(gca, 'Ylim', [Y_min Y_max]);
  
  hold off
  
end



function signal_flt = filter_eeg(signal, flt)
% Anton Tokariev (University of Helsinki, BABA Center, Finland)
%
% Input:
% signal: [samples x ch]
%    flt: cell array with filter objects (top - HPF, bottom - LPF)
%    Fs: sampling rate
%
% Output:
% signal_flt: cell array{1, N of fr.bands} of filtered data [samples x ch]
 
 % Signal Length
   L = size(signal, 1);

 % Flipped signal     
   signal_ud = flipud(signal);
   
 % Add pieces to signal (to fix edge effects)
   signal_big = [signal_ud; signal; signal_ud];                                                 
  
 % Number of bandpass filters (= fr.bands); columns = bandpass filters  
   N_flt = size(flt, 2); 
 
 % Init cell array for bandpass filtered signals   
   signal_flt{1, N_flt} = []; 
   
% Filter signals (bandpass filter = HPF + LPF)  
  for fr = 1:N_flt
      
      buf = [];                                                            %#ok<NASGU>
      buf = filtfilt(flt{1, fr}, signal_big); % HPF/cutoff = 0.85xFc 
      buf = filtfilt(flt{2, fr}, buf);        % LPF/cutoff = 1.15xFc
      
      signal_flt{1, fr} = buf(L+1:L+L, :);    % cut signal_big >> orig.sig.
      
  end
    
end



function signal_src = get_parcel_signals(signal, InverseOperator, CollapseOperator, MyAtlas)
 % Anton Tokariev (University of Helsinki, BABA Center, Finland)
 %
 % see also Tokariev et al., 2019, Cerebral Cortex
 %
 % Input:
 % signal: cell array {1 x N freq.} of filtered EEG data [samples x ch]
 % InverseOperator: Inverse solution for infant head with N = 19 trodes
 % CollapseOperator: weights for src signals within 'host' parcels
 % MyAtlas: assignment of src/verticies to parcels (in MyAtlas.Parcels)
 %
 % Output:
 % signal_src: cell array {1 x N freq.} of filtered parcel signals [samples x ch] 
 
   N_fr = size(signal, 2);       % number of frequencies
   
   Np = length(MyAtlas.Parcels); % number of parcels
   
   L = size(signal{1, 1}, 1);    % EEG length 
   
   signal_src{1, N_fr} = [];     % init output array
   
   CollapseOperator = repmat(CollapseOperator, 1, L);
   
   for k = 1:N_fr
      
     % source signals  
       src_buf = (InverseOperator * signal{1, k}').* CollapseOperator;     % [sources x samples]
       
     % parcel/cortical signals  
       parcel_buf = zeros(Np, L);
  
       for j = 1:Np
           parcel_buf(j, :) = mean(src_buf(MyAtlas.Parcels{j, 1}, :));     % [parcels x samples] 
       end
        
       signal_src{1, k} = parcel_buf';                                     % [samples x parcels] 
       
   end

end



function wpli = get_wPLI(A)
% This function computes debiased weighted phase lag index (dwPLI)
%
% It was adapted from the original implementation in Fieldtrip toolbox:
% www.fieldtriptoolbox.org/reference/ft_connectivity_wpli
% For more details read: Vinck et al., 2011, Neuroimage 

% Input:     A - band-pass filtered signals [samples x ch]
% Output: wpli - PPC (dwPLI) connectivity matrix [ch x ch]

  L = size(A, 1); % length
 ch = size(A, 2); % ch
 
 % Compute (debiased) wPLI
    wpli = zeros(ch, ch);
   
       for ChA = 1:ch
       
           X = A(:, ChA);     % signal channelA
        
           for ChB = ChA+1:ch
           
                Y = A(:, ChB); % signal channelB
           
                   Pxy = [];                                               %#ok<NASGU>  
                 % cross-spectral density
                   Pxy = cpsd(X, Y, L, 1, [], []);  
                 % compute wpli
                   Pxy = imag(Pxy); % make everything imaginary 
                   
                   outsum = [];                                            %#ok<NASGU>
                   outsum = nansum(Pxy, 1); % compute the sum
                   
                   outsumW = [];                                           %#ok<NASGU>
                   outsumW = nansum(abs(Pxy), 1); % normalization of WPLI 
                   
                  % Debiased version 
                    outssq = [];                                           %#ok<NASGU>
                    outssq = nansum(Pxy .^ 2, 1);
                    
                    wpli(ChA, ChB) = (outsum .^ 2 - outssq) ./ (outsumW .^ 2 - outssq); % do the pairwise    
                   
           end % ChB
           
       end % ChA
 
end % end



function visualize_head_model()

% This function plots:
% 1. scalp surface
% 2. EEG electrodes (red)
% 3. smoothed cortex
% 4. parcel centroids: Front./orange, Cent./purple, Occ./black, Temp./green

 
% load data
  load('Head Model\Atlas.mat');                                            %#ok<LOAD>
  load('Head Model\eeg_cap_19ch.mat');                                     %#ok<LOAD>
  load('Head Model\surfaces\cortex_smoothed.mat');                         %#ok<LOAD>
  load('Head Model\surfaces\scalp.mat');                                   %#ok<LOAD>
  
  Np = size(Atlas.Centroids, 1); % N of parcels
   
  figure;
   hold on
    set(gcf, 'Color', 'w');
    
   % plot surfaces 
     patch('Vertices', flat_cx.Vertices, 'Faces', flat_cx.Faces, 'FaceColor', [0.80 0.80 0.80], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
     patch('Vertices', scalp.Vertices,   'Faces', scalp.Faces,   'FaceColor', [0.93 0.84 0.84], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
   % colour coded parcel centroids
     scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 3), 40, 'ok', 'filled');
     scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 3), 40, 'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerEdgeColor', [0.47 0.67 0.19]);
     scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 3), 40, 'MarkerFaceColor', [0.75 0.00 0.75], 'MarkerEdgeColor', [0.75 0.00 0.75]);
     scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 3), 40, 'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', [0.93 0.69 0.13]);
  
   % centroid numbers  
     for j = 1:Np
         text(1.1*Atlas.Centroids(j, 1), 1.1*Atlas.Centroids(j, 2), 1.1*Atlas.Centroids(j, 3), num2str(j), 'FontSize', 7, 'Color', 'b', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
     end
   
   % eeg electrodes  
     scatter3(cap19(:, 1), cap19(:, 2), cap19(:, 3), 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); %#ok<USENS>
    
    axis off
    
    xlim([-0.07 0.07]);
    ylim([-0.06 0.06]);
    zlim([-0.04 0.09]);
    
    view([35 15]);
    
    title('Head model');
    
   hold off  

end


function show_3D_network(Matrix, Atlas, ttl)

% Input: Matrix - (binary) connectivity matrix [ch x ch]
%        Atlas.Centroids - Cartesian coordinate of parcel centroids (xyz)
%        Atlas.Areas - Parcel assingments to Frontal (F), Central (C),
%                      Occipital (O), and Temporal (T) brain areas
%        ttl - title (string)

 ch = size(Matrix, 1);
 
 figure;
  hold on
   set(gcf, 'Color', 'w');

 % plot network connections 
   for ChA = 1:ch
       for ChB = ChA+1:ch
           if Matrix(ChA, ChB) == 1
              line([Atlas.Centroids(ChA, 1) Atlas.Centroids(ChB, 1)], [Atlas.Centroids(ChA, 2) Atlas.Centroids(ChB, 2)], [Atlas.Centroids(ChA, 3) Atlas.Centroids(ChB, 3)], 'Color', 'r', 'Linewidth', 1); 
           end
       end
   end
  
 % plot colour coded parcel centroids
   scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'O', 3), 40, 'ok', 'filled');
   scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'T', 3), 40, 'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerEdgeColor', [0.47 0.67 0.19]);
   scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'C', 3), 40, 'MarkerFaceColor', [0.75 0.00 0.75], 'MarkerEdgeColor', [0.75 0.00 0.75]);
   scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'F', 3), 40, 'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', [0.93 0.69 0.13]);
 
 % plot smoothed cx  
   load('Head Model\surfaces\cortex_smoothed.mat'); %#ok<LOAD>
   patch('Vertices', flat_cx.Vertices, 'Faces', flat_cx.Faces, 'FaceColor', [0.80 0.80 0.80], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

   xlim([-0.06 0.06]);
   ylim([-0.05 0.05]);
   zlim([0 0.08]);
   
   axis off
    
   view([0 90]);
    
   title(ttl);
    
  hold off 
   
end


 
