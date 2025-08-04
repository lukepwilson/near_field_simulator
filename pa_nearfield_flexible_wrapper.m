% compute_chrank_heatmap_with_tx_parallel.m
% STATIC HEATMAP OF CHANNEL RANK vs. VERTICAL & HORIZONTAL OFFSETS
% WITH TX ANTENNA LOCATIONS AND REGION CONTOURS OVERLAID
% Parallelized version using parfor
% Luke Wilson
clear; clc; close all;

%% USER PARAMETERS
fc = 29e9; % center frequency (Hz)
N_pa_tx = 3; % # TX subarrays
N_pa_rx = 3; % # RX subarrays
spread_tx = 10; % TX per subarray spread (m)
spread_rx = 1; % RX per subarray spread (m)

% sweep ranges
vertical_offset_vals = linspace(0.01, 30, 50); % RX–TX separation (m)
horizontal_offset_vals = linspace(-80, 80, 50); % lateral offset (m)

% beam & null angles (broadside + 90° null)
tx_beam_angles = zeros(1, N_pa_tx);
tx_beam_angles(1) = 40;
tx_beam_angles(2) = 0;
tx_beam_angles(3) = -40;


rx_beam_angles = zeros(1, N_pa_rx);
rx_beam_angles(1) = 40;
rx_beam_angles(2) = 0;
rx_beam_angles(3) = -40;


tx_null_angles = 90 * ones(1, N_pa_tx);
rx_null_angles = 90 * ones(1, N_pa_rx);

%% SETUP PARALLEL POOL
% Check if parallel pool exists, create if needed
if isempty(gcp('nocreate'))
    parpool(); % Uses default number of workers
end

%% PREALLOCATE AND COMPUTE RANK MAP
M = numel(vertical_offset_vals);
N = numel(horizontal_offset_vals);
rank_map = nan(M, N);
total = M * N;

% Create linear indices for progress tracking
[I, J] = meshgrid(1:M, 1:N);
i_vals = I(:);
j_vals = J(:);

% Preallocate results vector for parfor
rank_results = zeros(total, 1);

tStart = tic;

% Parallel computation with linear indexing
parfor idx = 1:total
    i = i_vals(idx);
    j = j_vals(idx);
    
    dtrx = vertical_offset_vals(i);
    offset = horizontal_offset_vals(j);
    
    [~, ~, ch_rank] = pa_nearfield_flexible_core( ...
        fc, dtrx, ...
        N_pa_rx, N_pa_tx, ...
        spread_rx, spread_tx, ...
        offset, ...
        tx_beam_angles, tx_null_angles, ...
        rx_beam_angles, rx_null_angles);
    
    rank_results(idx) = ch_rank;  % Store in linear array
    
    % Progress reporting (note: this won't be ordered due to parallelization)
    if mod(idx, 100) == 0 || idx == total
        fprintf('Completed %d/%d iterations\n', idx, total);
    end
end

% Reshape results back to 2D matrix
rank_map = reshape(rank_results, M, N).';

elapsed = toc(tStart);
fprintf('Heavy computation took %.2f seconds.\n', elapsed);

%% ALTERNATIVE: Row-wise parfor approach (SIMPLER OPTION)
% If the above doesn't work, uncomment this section and comment out the above parfor
% This approach parallelizes over rows instead of all elements

% tStart = tic;
% parfor i = 1:M
%     dtrx = vertical_offset_vals(i);
%     temp_rank = zeros(1, N);
%     
%     for j = 1:N
%         offset = horizontal_offset_vals(j);
%         [~, ~, ch_rank] = pa_nearfield_flexible_core( ...
%             fc, dtrx, ...
%             N_pa_rx, N_pa_tx, ...
%             spread_rx, spread_tx, ...
%             offset, ...
%             tx_beam_angles, tx_null_angles, ...
%             rx_beam_angles, rx_null_angles);
%         temp_rank(j) = ch_rank;
%     end
%     
%     rank_map(i, :) = temp_rank;
%     fprintf('Completed row %d/%d\n', i, M);
% end
% elapsed = toc(tStart);
% fprintf('Heavy computation took %.2f seconds.\n', elapsed);

%% GET TX ELEMENT LOCATIONS
[~, tx_loc, ~] = pa_nearfield_flexible_core( ...
    fc, vertical_offset_vals(1), ...
    N_pa_rx, N_pa_tx, ...
    spread_rx, spread_tx, ...
    0, ... % no horizontal offset
    tx_beam_angles, tx_null_angles, ...
    rx_beam_angles, rx_null_angles);

%% DEFINE REGION BOUNDARIES
d_tx = (N_pa_tx - 1) * spread_tx;
d_rx = (N_pa_rx - 1) * spread_rx;
lambda = 3e8 / fc;
region_1 = d_tx * d_rx / lambda;
region_2 = 2 * region_1;

%% CREATE GRID FOR CONTOURS
[H_grid, V_grid] = meshgrid(horizontal_offset_vals, vertical_offset_vals);
radial_dist = sqrt(H_grid.^2 + V_grid.^2);

%% PLOT HEATMAP
figure;
imagesc(horizontal_offset_vals, vertical_offset_vals, rank_map);
set(gca, 'YDir', 'normal');
colormap(jet);
cb = colorbar;
clim([1, min(N_pa_tx, N_pa_rx)]);
xlabel('Horizontal Offset (m)');
ylabel('Vertical RX–TX Separation (m)');
title('Channel Rank Heatmap with Regions');
hold on;

% Overlay region boundaries
contour(H_grid, V_grid, radial_dist, [region_1 region_1], ...
    'LineColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Region 1');
contour(H_grid, V_grid, radial_dist, [region_2 region_2], ...
    'LineColor', 'b', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Region 2');

% Overlay TX antenna x-positions at y = 0
plot(tx_loc, zeros(size(tx_loc)), 'kp', ...
    'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', 'TX elements');

legend('Location', 'northeastoutside');
hold off;

%% OPTIONAL: Save results for future use
% save('channel_rank_heatmap_results.mat', 'rank_map', 'vertical_offset_vals', ...
%      'horizontal_offset_vals', 'tx_loc', 'region_1', 'region_2');