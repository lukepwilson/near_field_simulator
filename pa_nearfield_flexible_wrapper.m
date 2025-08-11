% % compute_chrank_heatmap_with_tx_parallel.m
% % STATIC HEATMAP OF CHANNEL RANK vs. VERTICAL & HORIZONTAL OFFSETS
% % WITH TX ANTENNA LOCATIONS AND REGION CONTOURS OVERLAID
% % Parallelized version using parfor
% % Luke Wilson
% clear; clc; close all;
% 
% %% USER PARAMETERS
% fc = 29e9; % center frequency (Hz)
% N_pa_tx = 2; % # TX subarrays
% N_pa_rx = 2; % # RX subarrays
% spread_tx = 1; % TX per subarray spread (m)
% spread_rx = 1; % RX per subarray spread (m)
% vertical_dist = 50;
% horizontal_dist = 50;
% number_samples_per_line = 100;
% 
% % sweep ranges
% vertical_offset_vals = linspace(0.01, vertical_dist, number_samples_per_line); % RX–TX separation (m)
% horizontal_offset_vals = linspace(-horizontal_dist, horizontal_dist, number_samples_per_line); % lateral offset (m)
% 
% % beam & null angles (broadside + 90° null)
% tx_beam_angles = zeros(1, N_pa_tx);
% tx_beam_angles(1) = 20;
% tx_beam_angles(2) = 0;
% tx_beam_angles(3) = -20;
% 
% 
% rx_beam_angles = zeros(1, N_pa_rx);
% rx_beam_angles(1) = 20;
% rx_beam_angles(2) = 0;
% rx_beam_angles(3) = -20;
% 
% 
% tx_null_angles = 90 * ones(1, N_pa_tx);
% rx_null_angles = 90 * ones(1, N_pa_rx);
% 
% %% SETUP PARALLEL POOL
% % Check if parallel pool exists, create if needed
% if isempty(gcp('nocreate'))
%     parpool(); % Uses default number of workers
% end
% 
% %% PREALLOCATE AND COMPUTE RANK MAP
% M = numel(vertical_offset_vals);
% N = numel(horizontal_offset_vals);
% rank_map = nan(M, N);
% total = M * N;
% 
% % Create linear indices for progress tracking
% [I, J] = meshgrid(1:M, 1:N);
% i_vals = I(:);
% j_vals = J(:);
% 
% % Preallocate results vector for parfor
% rank_results = zeros(total, 1);
% 
% tStart = tic;
% 
% % Parallel computation with linear indexing
% parfor idx = 1:total
%     i = i_vals(idx);
%     j = j_vals(idx);
% 
%     dtrx = vertical_offset_vals(i);
%     offset = horizontal_offset_vals(j);
% 
%     [~, ~, ch_rank] = pa_nearfield_flexible_core( ...
%         fc, dtrx, ...
%         N_pa_rx, N_pa_tx, ...
%         spread_rx, spread_tx, ...
%         offset, ...
%         tx_beam_angles, tx_null_angles, ...
%         rx_beam_angles, rx_null_angles);
% 
%     rank_results(idx) = ch_rank;  % Store in linear array
% 
%     % Progress reporting (note: this won't be ordered due to parallelization)
%     if mod(idx, 100) == 0 || idx == total
%         fprintf('Completed %d/%d iterations\n', idx, total);
%     end
% end
% 
% % Reshape results back to 2D matrix
% rank_map = reshape(rank_results, M, N).';
% 
% elapsed = toc(tStart);
% fprintf('Heavy computation took %.2f seconds.\n', elapsed);
% 
% %% ALTERNATIVE: Row-wise parfor approach (SIMPLER OPTION)
% % If the above doesn't work, uncomment this section and comment out the above parfor
% % This approach parallelizes over rows instead of all elements
% 
% % tStart = tic;
% % parfor i = 1:M
% %     dtrx = vertical_offset_vals(i);
% %     temp_rank = zeros(1, N);
% %     
% %     for j = 1:N
% %         offset = horizontal_offset_vals(j);
% %         [~, ~, ch_rank] = pa_nearfield_flexible_core( ...
% %             fc, dtrx, ...
% %             N_pa_rx, N_pa_tx, ...
% %             spread_rx, spread_tx, ...
% %             offset, ...
% %             tx_beam_angles, tx_null_angles, ...
% %             rx_beam_angles, rx_null_angles);
% %         temp_rank(j) = ch_rank;
% %     end
% %     
% %     rank_map(i, :) = temp_rank;
% %     fprintf('Completed row %d/%d\n', i, M);
% % end
% % elapsed = toc(tStart);
% % fprintf('Heavy computation took %.2f seconds.\n', elapsed);
% 
% %% GET TX ELEMENT LOCATIONS
% [~, tx_loc, ~] = pa_nearfield_flexible_core( ...
%     fc, vertical_offset_vals(1), ...
%     N_pa_rx, N_pa_tx, ...
%     spread_rx, spread_tx, ...
%     0, ... % no horizontal offset
%     tx_beam_angles, tx_null_angles, ...
%     rx_beam_angles, rx_null_angles);
% 
% %% DEFINE REGION BOUNDARIES
% d_tx = (N_pa_tx - 1) * spread_tx;
% d_rx = (N_pa_rx - 1) * spread_rx;
% lambda = 3e8 / fc;
% region_1 = d_tx * d_rx / lambda;
% region_2 = 2 * region_1;
% 
% %% CREATE GRID FOR CONTOURS
% [H_grid, V_grid] = meshgrid(horizontal_offset_vals, vertical_offset_vals);
% radial_dist = sqrt(H_grid.^2 + V_grid.^2);
% 
% %% PLOT HEATMAP
% figure;
% imagesc(horizontal_offset_vals, vertical_offset_vals, rank_map);
% set(gca, 'YDir', 'normal');
% colormap(jet);
% cb = colorbar;
% clim([1, min(N_pa_tx, N_pa_rx)]);
% xlabel('Horizontal Offset (m)');
% ylabel('Vertical RX–TX Separation (m)');
% title('Channel Rank Heatmap with Regions');
% hold on;
% 
% % Overlay region boundaries
% contour(H_grid, V_grid, radial_dist, [region_1 region_1], ...
%     'LineColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Region 1');
% contour(H_grid, V_grid, radial_dist, [region_2 region_2], ...
%     'LineColor', 'b', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Region 2');
% 
% % Overlay TX antenna x-positions at y = 0
% plot(tx_loc, zeros(size(tx_loc)), 'kp', ...
%     'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', 'TX elements');
% 
% legend('Location', 'northeastoutside');
% hold off;
% 
% %% OPTIONAL: Save results for future use
% % save('channel_rank_heatmap_results.mat', 'rank_map', 'vertical_offset_vals', ...
% %      'horizontal_offset_vals', 'tx_loc', 'region_1', 'region_2');

% compute_chrank_heatmap_with_tx_parallel.m
% STATIC HEATMAP (FUSION): Hue encodes CHANNEL RANK, Brightness encodes σ1
% WITH TX ANTENNA LOCATIONS AND REGION CONTOURS OVERLAID
% Parallelized version using parfor
% Luke Wilson (fusion rendering added)

%% ====================== USER PARAMETERS ======================
fc   = 29e9;        % center frequency (Hz)
N_pa_tx = 2;        % # TX subarrays
N_pa_rx = 2;        % # RX subarrays
spread_tx = 1;      % TX per subarray spread (m)
spread_rx = 1;      % RX per subarray spread (m)
vertical_dist   = 2000;
horizontal_dist = 50;
number_samples_per_line = 50;

% sweep ranges
vertical_offset_vals   = linspace(0.01, vertical_dist,   number_samples_per_line); % RX–TX separation (m)
horizontal_offset_vals = linspace(-horizontal_dist, horizontal_dist, number_samples_per_line); % lateral offset (m)

% beam angles (simple pattern: [20, 0, -20] trimmed/extended to N)
tx_beam_angles = [-10 10];
rx_beam_angles = [10 -10];


% null angles (broadside + 90° null)
tx_null_angles = 90 * ones(1, N_pa_tx);
rx_null_angles = 90 * ones(1, N_pa_rx);

%% ====================== PARALLEL POOL ========================
if isempty(gcp('nocreate'))
    parpool(); % Uses default number of workers
end

%% ====================== PREALLOCATE / INDICES ================
M = numel(vertical_offset_vals);
N = numel(horizontal_offset_vals);
total = M * N;

[I, J] = meshgrid(1:M, 1:N);
i_vals = I(:);
j_vals = J(:);

rank_results   = zeros(total, 1);
sigma1_results = zeros(total, 1);

tStart = tic;

%% ====================== PARFOR: CORE SWEEP ===================
parfor idx = 1:total
    i = i_vals(idx);
    j = j_vals(idx);

    dtrx  = vertical_offset_vals(i);
    offset = horizontal_offset_vals(j);

    [eig_vals, ~, ch_rank] = pa_nearfield_flexible_core( ...
        fc, dtrx, ...
        N_pa_rx, N_pa_tx, ...
        spread_rx, spread_tx, ...
        offset, ...
        tx_beam_angles, tx_null_angles, ...
        rx_beam_angles, rx_null_angles);

    % First singular value magnitude
    s1 = max(abs(eig_vals));     % assume eig_vals are singular values or related magnitudes
    sigma1_results(idx) = s1;
    rank_results(idx)   = ch_rank;

    % (Optional) sparse progress print
    if mod(idx, 100) == 0
        fprintf('Completed %d/%d iterations\n', idx, total);
    end
end

elapsed = toc(tStart);
fprintf('Heavy computation took %.2f seconds.\n', elapsed);

% Reshape to grids (match imagesc/image orientation expectations)
rank_map   = reshape(rank_results,   M, N).';
sigma1_map = reshape(sigma1_results, M, N).';

%% ====================== TX LOCATIONS (y=0 overlay) ===========
[~, tx_loc, ~] = pa_nearfield_flexible_core( ...
    fc, vertical_offset_vals(1), ...
    N_pa_rx, N_pa_tx, ...
    spread_rx, spread_tx, ...
    0, ... % no horizontal offset
    tx_beam_angles, tx_null_angles, ...
    rx_beam_angles, rx_null_angles);

%% ====================== REGION BOUNDARIES ====================
d_tx  = (N_pa_tx - 1) * spread_tx;
d_rx  = (N_pa_rx - 1) * spread_rx;
lambda = 3e8 / fc;
region_1 = d_tx * d_rx / lambda;
region_2 = 2 * region_1;

% Grid for contours
[H_grid, V_grid] = meshgrid(horizontal_offset_vals, vertical_offset_vals);
radial_dist = sqrt(H_grid.^2 + V_grid.^2);

%% ====================== FUSION: H=RANK, V=σ1 =================
Rmax = min(N_pa_tx, N_pa_rx);
rank_clipped = max(1, min(Rmax, round(rank_map)));

% σ1 in dB -> normalize to [0,1] for brightness
sigma1_db = 20*log10(sigma1_map + eps);

% Robust range: percentile clamp to avoid outliers washing out image
p_lo = 5; p_hi = 95;
db_lo = prctile(sigma1_db(~isinf(sigma1_db(:)) & ~isnan(sigma1_db(:))), p_lo);
db_hi = prctile(sigma1_db(~isinf(sigma1_db(:)) & ~isnan(sigma1_db(:))), p_hi);
if isempty(db_lo) || isempty(db_hi) || db_lo == db_hi
    db_lo = min(sigma1_db(:));
    db_hi = max(sigma1_db(:));
end
V = (sigma1_db - db_lo) ./ max(db_hi - db_lo, eps);
V = min(max(V,0),1);

% Hue mapping: 1..Rmax → blue→red (HSV hue 0.67→0.0)
hue_start = 0.67; hue_end = 0.0;
H = hue_start + (hue_end - hue_start) * (rank_clipped - 1) / max(Rmax - 1, 1);

% Full saturation for vivid hue; fusion uses brightness to encode σ1
S = ones(size(H));

HSV = cat(3, H, S, V);
RGB = hsv2rgb(HSV);

%% ====================== PLOT: FUSED IMAGE + OVERLAYS =========
fig = figure('Name','Rank (Hue) + \sigma_1 (Brightness)','NumberTitle','off');
axMain = axes('Parent',fig);
image(axMain, horizontal_offset_vals, vertical_offset_vals, RGB);
set(axMain, 'YDir', 'normal');
xlabel(axMain, 'Horizontal Offset (m)');
ylabel(axMain, 'Vertical RX–TX Separation (m)');
title(axMain, 'Rank (Hue) + First Singular Value (Brightness)');
hold(axMain, 'on');

% Region contours - with error handling
hR1 = contour(axMain, H_grid, V_grid, radial_dist, [region_1 region_1], ...
    'LineColor','k','LineWidth',1.5,'DisplayName','Region 1');
hR2 = contour(axMain, H_grid, V_grid, radial_dist, [region_2 region_2], ...
    'LineColor','b','LineStyle','--','LineWidth',1.5,'DisplayName','Region 2');

% TX elements at y=0
hTX = plot(axMain, tx_loc, zeros(size(tx_loc)), 'kp', ...
    'MarkerFaceColor','w', 'MarkerSize',8, 'DisplayName','TX elements');

% Discrete legend entries for ranks 1..Rmax
rank_colors = hsv2rgb([linspace(hue_start, hue_end, Rmax).'  ones(Rmax,1)  ones(Rmax,1)]);
hRank = gobjects(Rmax,1);
for r = 1:Rmax
    hRank(r) = plot(axMain, nan, nan, 's', 'MarkerSize',10, ...
        'MarkerFaceColor', rank_colors(r,:), 'MarkerEdgeColor','k', ...
        'DisplayName', sprintf('Rank %d', r));
end

% Build legend with available handles
legendHandles = [];
legendLabels = {};

% Check if contour handles are valid and non-empty
if ~isempty(hR1) && ishandle(hR1(1))
    legendHandles(end+1) = hR1(1);
    legendLabels{end+1} = 'Region 1';
else
    % Create dummy handle if contour is out of range
    hR1_dummy = plot(axMain, nan, nan, 'k-', 'LineWidth', 1.5);
    legendHandles(end+1) = hR1_dummy;
    legendLabels{end+1} = sprintf('Region 1 (r=%.1fm)', region_1);
end

if ~isempty(hR2) && ishandle(hR2(1))
    legendHandles(end+1) = hR2(1);
    legendLabels{end+1} = 'Region 2';
else
    % Create dummy handle if contour is out of range
    hR2_dummy = plot(axMain, nan, nan, 'b--', 'LineWidth', 1.5);
    legendHandles(end+1) = hR2_dummy;
    legendLabels{end+1} = sprintf('Region 2 (r=%.1fm)', region_2);
end

% Add TX handle
legendHandles(end+1) = hTX;
legendLabels{end+1} = 'TX elements';

% Add rank handles
for r = 1:Rmax
    legendHandles(end+1) = hRank(r);
    legendLabels{end+1} = sprintf('Rank %d', r);
end

% Create legend with collected handles
legend(legendHandles, legendLabels, 'Location','northeastoutside');

% Tight layout
axis(axMain, 'tight');

%% ====================== OPTIONAL: SAVE RESULTS ================
% save('channel_rank_fusion_results.mat', 'rank_map','sigma1_map', ...
%      'vertical_offset_vals','horizontal_offset_vals', ...
%      'tx_loc','region_1','region_2','db_lo','db_hi');

%% ====================== LOCAL FUNCTIONS ======================
function angs = default_beam_angles(N)
% Returns a simple symmetric set of beam angles sized to N
    base = [20, 0, -20];
    if N <= numel(base)
        angs = base(1:N);
    else
        % Extend by padding zeros if more than 3; adjust as desired
        angs = zeros(1,N);
        k = min(numel(base), N);
        angs(1:k) = base(1:k);
    end
end