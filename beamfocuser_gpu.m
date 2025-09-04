% compute_chrank_heatmap_fast.m
% Fast batched GPU version (fallback to CPU if no GPU)
% Luke Wilson

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%%%  ARRAY GEOMETRY
%%%%%%%%%%%%%%%%%%%
fc = 29e9;              % center frequency (Hz)
N_pa_tx = 100;            % number of TX subarrays
spread_tx = 1000;         % TX per subarray spread (m)
n_ant = 36;             % number of antennas per subarray
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%%%  PLOT CHARACTERISTICS
%%%%%%%%%%%%%%%%%%%
max_vertical_dist   = 100000e3;     % maximum vertical viewing distance (m)
min_vertical_dist = 0.01;
max_horizontal_dist = 200e3;     % maximum horizontal viewing distance (m)
log_scale_axis = 1;             % pick whether you want the y-axis to be in linear or log scale
ant_to_plot = 1:N_pa_tx;         % index of phased arrays to plot.
number_of_x_points = 200;
number_of_y_points = 200;

min_svd = -80;
max_svd = 0;

min_capacity = 30;
max_capacity = 1000;

%%%%%%%%%%%%%%%%%%%
%%%  Beam Focusing
%%%%%%%%%%%%%%%%%%%
beam_angle    = 0;      % deg
beam_distance = 2000e3;    % m
phase_alignment = 1;

%%%%%%%%%%%%%%%%%%%
%%%  THROUGHPUT / NOISE MODEL
%%%%%%%%%%%%%%%%%%%
B_Hz          = 100e6;      % bandwidth
Tx_power_W    = 1;          % total TX power (watts)
T_kelvin      = 290;        % noise temp
noiseFig_dB   = 5;          % receiver noise figure (dB)
k_B           = 1.380649e-23; % Boltzmann
N_r           = N_pa_tx;      % #Rx chains (>= #streams you want)
rx_spacing    = 3e8/(fc*2);        % UE antenna spacing (e.g., λ/2)
use_MRT_across_subarrays = false; % if true, treat as 1 stream (beamforming), not multiplexing
%%%%%%%%%%%%%%%%%%%


% sweep ranges (usually don't need to change)
horizontal_offset_vals  = linspace(-max_horizontal_dist, max_horizontal_dist, number_of_x_points); % X axis (linear)

if log_scale_axis
    vertical_offset_vals    = logspace(log10(min_vertical_dist), log10(max_vertical_dist), number_of_y_points);      % Y axis (log)
else
    vertical_offset_vals    = linspace(min_vertical_dist, max_vertical_dist, number_of_y_points);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXECUTION MODE: GPU vs CPU
useGPU = parallel.gpu.GPUDevice.isAvailable;
if useGPU, gpuDevice([]); end

lambda = 3e8 / fc;
d_int  = lambda/2;

%% ---------- Build TX geometry & w_tx ONCE (CPU), then cast to GPU if needed ----------
% TX element x-positions (centered)
tx_phased_arrays_loc = zeros(N_pa_tx, n_ant);

for i = 1:N_pa_tx
    tx_phased_arrays_loc(i,:) = (0:n_ant-1) * d_int;   
end

tx_spread = tx_phased_arrays_loc + (0:N_pa_tx-1)' * spread_tx; % [N_pa_tx x n_ant]
tx_flat   = reshape(tx_spread.',1,[]);                                     % [1 x (N_pa_tx*n_ant)]
tx_flat   = tx_flat - mean(tx_flat);
num_tx    = numel(tx_flat);

% Per-subarray pointing angles -> digital weights
poi_x = beam_distance * sind(beam_angle);
poi_y = beam_distance * cosd(beam_angle);
avg_dist = zeros(1,N_pa_tx);
avg_angle = zeros(1, N_pa_tx);

% if phase_alignment
%     distances_for_phases = zeros(N_pa_tx,n_ant); 
%     phases_for_weights = zeros(N_pa_tx,n_ant);
% end

for i = 1:N_pa_tx
    i0 = (i-1)*n_ant + 1; 
    i1 = i*n_ant;
    x1 = tx_flat(i0);    x2 = tx_flat(i1);
    y1 = 0;              y2 = 0;
    ang1 = rad2deg(atan2(poi_y - y1, poi_x - x1));
    ang2 = rad2deg(atan2(poi_y - y2, poi_x - x2));
    dist1 = sqrt( (poi_y - y1).^2 + (poi_x - x1).^2 );
    dist2 = sqrt( (poi_y - y2).^2 + (poi_x - x2).^2 );
    avg_angle(i) = mean([ang1, ang2]);
    avg_dist(i) = mean([dist1,dist2]);
    % if phase_alignment
    %     for ii = 1:n_ant
    %         dx = poi_x - tx_flat(n_ant*(i-1)+ii);
    %         dy = poi_y - 0;
    %         distances_for_phases(i,ii) = sqrt( dx.^2 + dy.^2 );
    %         phases_for_weights(i,ii) = exp(1j * 2*pi .*  distances_for_phases(i,ii)./ single(lambda));
    %     end
    % end

end
for ii = 1:numel(avg_angle)
    if     avg_angle(ii) > 0,  avg_angle(ii) =  90 - avg_angle(ii);
    elseif avg_angle(ii) < 0,  avg_angle(ii) = -90 - avg_angle(ii);
    else,                      avg_angle(ii) =   0;
    end
end




w_tx = zeros(N_pa_tx, num_tx, 'single');   % single precision is faster/stable on GPU
psuedo_angles = zeros(1,N_pa_tx);
for i = ant_to_plot
    [w_row,psuedo_angles(i)] = get_singlefocal_weights(avg_angle(i),avg_dist(i),n_ant,fc,1,0); % [1 x n_ant] (double)
    i0 = (i-1)*n_ant + 1; i1 = i*n_ant;
    % if phase_alignment
    %     w_tx(i, i0:i1) = single(w_row .*  phases_for_weights(i,:).');        % cast to single4
    % else
        w_tx(i, i0:i1) = single(w_row);
    % end

    %% ----- Per-subarray SVD at the focal point -----
    
% Build per-element channel from TX element -> focal point
DXf = (poi_x - tx_flat(:));              % [num_tx x 1]
DYf = (poi_y - 0);                       % scalar (all TX at y=0 here)
Rf  = sqrt(DXf.^2 + DYf.^2);             % [num_tx x 1]
h_elem_focus = (1 ./ sqrt(4*pi*Rf)) .* exp(-1j * 2*pi * Rf / lambda);  % [num_tx x 1]

% Effective complex field at the focus per subarray (scalar each)
% (w_tx is [N_pa_tx x num_tx] with each row nonzero only on its own element block)
E_focus = w_tx * h_elem_focus;           % [N_pa_tx x 1] complex

reference_angle = angle(E_focus(1));
current_angle = angle( E_focus(i));

if i ~= 0
    w_tx(i,:) = w_tx(i,:) .* exp( 1j * ( reference_angle - current_angle ) );
end


% "SVD value" per subarray at the focal point = |scalar|
sv_per_array = abs(E_focus);             % [N_pa_tx x 1]
sv_per_array_dB = 20*log10(max(sv_per_array, realmin));

% Joint (all-subarrays combined) first singular value at the focal point
% (vector's 2-norm equals first singular value)
sv_joint = norm(E_focus);
sv_joint_dB = 20*log10(max(sv_joint, realmin));

end

if phase_alignment
    E_focus = w_tx * h_elem_focus;
end



% function C_bpsHz = su_mimo_capacity_at_point(x0, y0, tx_flat, w_tx, N_pa_tx, N_r, rx_spacing, lambda, B_Hz, Tx_power_W, T_kelvin, noiseFig_dB, k_B, useGPU,use_MRT_across_subarrays)
%     % Build Rx antenna coords (ULA centered at x0)
%     xr = x0 + ( (0:N_r-1) - (N_r-1)/2 ) * rx_spacing;
%     yr = y0 * ones(1, N_r);
% 
%     % Per-element channels to each Rx antenna: H_elem_full [num_tx x N_r]
%     num_tx = numel(tx_flat);
%     XE = reshape(tx_flat(:), [num_tx,1]);
%     YE = zeros(num_tx,1);
% 
%     XR = reshape(xr(:).', [1, N_r]);
%     YR = reshape(yr(:).', [1, N_r]);
% 
%     DX = XR - XE;                            % [num_tx x N_r]
%     DY = YR - YE;
%     R  = sqrt(DX.^2 + DY.^2);
%     H_elem_full = (1 ./ sqrt(4*pi*R)) .* exp(-1j * 2*pi .* R ./ lambda);   % [num_tx x N_r]
% 
%     % Combine within each subarray using your digital weights: W [N_pa_tx x num_tx]
%     % Effective channel from subarrays to Rx antennas: H_eff [N_r x N_pa_tx]
%     H_eff = (w_tx * H_elem_full).';   % (N_pa_tx x N_r)^T => [N_r x N_pa_tx]
% 
%     % Noise power with NF:
%     N0B_W = k_B * T_kelvin * B_Hz * 10^(noiseFig_dB/10);
%     rho   = Tx_power_W / N0B_W;       % total SNR (linear)
% 
%     if use_MRT_across_subarrays
%         % 1-stream beamforming upper bound (not multiplexing): equal power across subarrays not used.
%         g = vecnorm(H_eff,2,1);              % per-subarray gains if you were to combine? (not needed)
%         % For MRT across subarrays (single stream), effective gain is largest singular value of H_eff:
%         s = svd(H_eff);
%         C_bpsHz = log2( 1 + rho * (s(1)^2) );   % all power to strongest mode
%     else
%         % Equal power across N_s = N_pa_tx independent streams:
%         Ns = size(H_eff,2);
%         C_bpsHz = real(log2(det( eye(N_r) + (rho/Ns) * (H_eff*H_eff') )));
%     end
% end


%% ---------- GPU batched path (no loops) ----------
if useGPU
    tStart = tic;

    % Make GPU arrays in single precision
    Xtx = gpuArray(single(tx_flat));                 % [1 x num_tx]
    Xtx = reshape(Xtx, [num_tx 1 1]);                % [num_tx x 1 x 1]
    Ytx = gpuArray(single(zeros(size(Xtx), 'like', Xtx))); % all TX at y=0

    % Build RX grids on GPU (broadcast to pages)
    [Xrx_grid, Yrx_grid] = meshgrid( ...
        gpuArray(single(horizontal_offset_vals)), ...
        gpuArray(single(vertical_offset_vals)) );
    % Xrx_grid, Yrx_grid: [M x N] each; make them [1 x M x N] for broadcasting
    Xrx = reshape(Xrx_grid, [1 size(Xrx_grid,1) size(Xrx_grid,2)]); % [1 x M x N]
    Yrx = reshape(Yrx_grid, [1 size(Yrx_grid,1) size(Yrx_grid,2)]); % [1 x M x N]

    % Distances from every TX elem to every RX point: [num_tx x M x N]
    DX = Xrx - Xtx;      % implicit expansion
    DY = Yrx - Ytx;
    R  = sqrt(DX.^2 + DY.^2);

    amp   = 1 ./ sqrt(4*pi*R);
    phase = exp(-1j * 2*pi .* R ./ single(lambda));
    H_elem = amp .* phase;                       % [num_tx x M x N], complex single

    % Multiply by digital combiner across all pages:
    % reshape w_tx to [N_pa_tx x num_tx x 1 x 1], H_elem to [num_tx x 1 x M x N]
    W = gpuArray(single(w_tx));
    W = reshape(W, [N_pa_tx, num_tx, 1, 1]);
    H = reshape(H_elem, [num_tx, 1, size(H_elem,2), size(H_elem,3)]);

    % Batched matmul: result [N_pa_tx x 1 x M x N]
    H_dig = pagemtimes(W, H);

    % First singular value of a column vector == its 2-norm.
    % So s1 = sqrt(sum(|H_dig|.^2, over dim=1 (rows))) -> keep pages [1 x 1 x M x N]
    s1 = sqrt(sum(abs(H_dig).^2, 1));      % [1 x 1 x M x N]
    svd_map = gather(squeeze(s1));         % [M x N] on CPU

    fprintf('GPU batched computation: %.2f s\n', toc(tStart));

else
    %% ---------- CPU fallback (parfor), still uses norm() instead of svd ----------
    p = gcp('nocreate'); if isempty(p), parpool; end
    M = numel(vertical_offset_vals);
    N = numel(horizontal_offset_vals);
    total = M*N;

    tStart = tic;
    q = parallel.pool.DataQueue; afterEach(q, @progressListener);
    send(q, struct('init',true,'total',total,'tStart',tStart));

    svd_results = zeros(total, 1);
    parfor idx = 1:total
        iterStart = tic;
        [i, j] = ind2sub([M, N], idx);
        dtrx   = vertical_offset_vals(i);
        offset = horizontal_offset_vals(j);

        % Recompute the per-point element channel on CPU (double)
        DX = (offset - tx_flat);   % 1 x num_tx
        DY = (dtrx  - 0);          % scalar
        R  = sqrt(DX.^2 + DY.^2);
        H_elem = (1 ./ sqrt(4*pi*R)) .* exp(-1j * 2*pi .* R ./ lambda); % 1 x num_tx
        H_dig  = w_tx * H_elem.';   % [N_pa_tx x 1] (w_tx is single; cast promotes)
        svd_results(idx) = norm(H_dig);  % first singular value

        send(q, toc(iterStart));
    end
    svd_map = reshape(svd_results, [numel(vertical_offset_vals), numel(horizontal_offset_vals)]);
    fprintf('CPU parfor computation: %.2f s\n', toc(tStart));
end

% ===== Compute per-point SNR and throughput maps =====
% Noise power N0B with noise figure:
N0B_W   = k_B * T_kelvin * B_Hz * db2lin(noiseFig_dB);  % watts
rho_lin = Tx_power_W / N0B_W;                            % total SNR (linear)

% We modeled equal power across N_pa_tx "TX chains" -> rho/N_pa_tx
% Effective channel power gain per point is ||h||^2 = |s1|^2
if useGPU
    gain_map = gather(abs(s1).^2);            % [M x N], ||h||^2
else
    gain_map = abs(svd_map).^2;               % [M x N], ||h||^2
end

SNR_map_lin       = (rho_lin / N_pa_tx) .* gain_map;  % per-point post-beamforming SNR
C_map_bps_per_Hz  = log2(1 + SNR_map_lin);            % spectral efficiency
C_map_bps         = B_Hz .* C_map_bps_per_Hz;         % throughput (bits/s)


%% ---------- Plot ----------
svd_map_pos = max(real(svd_map), realmin);
svd_map_pos(~isfinite(svd_map_pos)) = realmin;
svd_map_db  = 20*log10(svd_map_pos);

figure('Color','w');
imagesc(horizontal_offset_vals, vertical_offset_vals, svd_map_db);

if log_scale_axis
    set(gca,'YDir','normal','YScale','log');
else
    set(gca,'YDir','normal');
end
% % set(gca,'YDir','normal','YScale','log');
% set(gca,'YDir','normal');
colormap(jet); cb = colorbar; cb.Label.String = 'First singular value (dB)';
clim([min_svd max_svd]);
xlabel('Horizontal Offset (m)'); ylabel('Vertical RX–TX Separation (m)');
title(sprintf('First Singular Value Heatmap (%s)', ternary(useGPU,'GPU','CPU')));

% Overlay TX x-locations at smallest Y (just for reference)
y_overlay = min(vertical_offset_vals);
hold on;
plot(tx_flat(:), y_overlay*ones(numel(tx_flat),1), 'kp', ...
    'MarkerFaceColor','w','MarkerSize',8,'DisplayName','TX elements');

% Overlay focal point
plot(poi_x, poi_y, 'ro', 'MarkerSize', 10, 'LineWidth', 2, ...
    'DisplayName','Focal Point');

legend('Location','northeastoutside');
hold off;


%% ---------- Capacity (fast) ----------
fprintf('\n=== Starting Capacity Computation (fast) ===\n');
tCapStart = tic;

% Grid vectors (reuse yours)
xv = linspace(-max_horizontal_dist, max_horizontal_dist, number_of_x_points);
if log_scale_axis
    yv = logspace(log10(min_vertical_dist), log10(max_vertical_dist), number_of_y_points);
else
    yv = linspace(min_vertical_dist, max_vertical_dist, number_of_y_points);
end
Nx = numel(xv); Ny = numel(yv); NT = Nx*Ny;

% Precompute RX antenna offsets (centered ULA around x0)
rx_offsets   = ((0:N_r-1) - (N_r-1)/2) * rx_spacing;   % 1 x N_r  (double)
rx_offsets_s = single(rx_offsets);

% Noise/SNR constants
N0B_W = single(k_B * T_kelvin * B_Hz * db2lin(noiseFig_dB));
rho   = single(Tx_power_W) / N0B_W;          % total SNR
Ns    = single(N_pa_tx);

% Output
C_map = zeros(Ny, Nx, 'single');

if useGPU
    % -------- GPU-assisted batched path WITH PROGRESS --------
    tx_flat_s = gpuArray(single(tx_flat(:)));            % [num_tx x 1]
    w_tx_s    = gpuArray(single(w_tx));                  % [N_pa_tx x num_tx]
    lambda_s  = single(lambda);

    % batching over grid points
    batch_size = 512;                                    % tune for VRAM
    [X0grid, Y0grid] = meshgrid(single(xv), single(yv)); % [Ny x Nx]
    X0v = reshape(X0grid.', [], 1);                      % NT x 1
    Y0v = reshape(Y0grid.', [], 1);                      % NT x 1

    total_batches = ceil(NT / batch_size);
    fprintf('[GPU Capacity] total points: %d | batch size: %d | batches: %d\n', NT, batch_size, total_batches);

    done_pts = 0;
    for b = 1:total_batches
        tBatch = tic;
        s = (b-1)*batch_size + 1;
        e = min(b*batch_size, NT);
        K = e - s + 1;

        X0_batch = X0v(s:e);    % [K x 1] single
        Y0_batch = Y0v(s:e);    % [K x 1] single

        % RX coords for all K points
        XR = gpuArray( reshape(X0_batch, [1 1 K]) ) + gpuArray( reshape(rx_offsets_s, [1 N_r 1]) ); % [1 x N_r x K]
        YR = gpuArray( reshape(Y0_batch, [1 1 K]) ) + gpuArray( zeros(1,N_r,1,'single','gpuArray') ); % [1 x N_r x K]

        % TX coords (all y=0)
        XE = reshape(tx_flat_s, [], 1, 1);    % [num_tx x 1 x 1]
        YE = gpuArray(zeros(size(XE), 'like', XE));

        % Distances -> element channel [num_tx x N_r x K]
        DX = XR - XE;
        DY = YR - YE;
        R  = sqrt(DX.^2 + DY.^2);
        H_elem = (1 ./ sqrt(4*pi*R)) .* exp(-1j * 2*pi .* R ./ lambda_s);

        % Reduce with digital weights to subarray dimension
        W = reshape(w_tx_s, [N_pa_tx, size(H_elem,1), 1]);
        H_eff_gpu = pagemtimes(W, H_elem);         % [N_pa_tx x N_r x K]
        H_eff_gpu = permute(H_eff_gpu, [2 1 3]);   % -> [N_r x N_pa_tx x K]

        % Gather small pages and do tiny SVD per page on CPU
        H_eff_cpu = gather(H_eff_gpu);             % (N_r x N_pa_tx x K), small
        C_batch = zeros(K,1,'single');

        % Small SVDs are super fast; still parallelize if you want
        parfor kk = 1:K
            svals = svd(double(H_eff_cpu(:,:,kk)), 'econ'); %#ok<PFBNS>
            sig2  = single(svals).^2;
            C_batch(kk) = sum(log2(1 + (rho/Ns) * sig2));  % b/s/Hz
        end

        % write back into C_map (we linearized with X fastest)
        C_map_tmp = reshape(C_batch, [Nx, 1]).';
        row_start = floor((s-1)/Nx)+1;
        col_start = mod((s-1),Nx)+1;

        idx = 1; r = row_start; c = col_start;
        while idx <= K
            take = min(Nx - c + 1, K - idx + 1);
            C_map(r, c:(c+take-1)) = C_map_tmp(1, idx:(idx+take-1));
            idx = idx + take;
            r = r + (c+take-1==Nx);
            c = 1;
        end

        % ---- progress print for this batch ----
        done_pts = done_pts + K;
        elapsed  = toc(tCapStart);
        batch_t  = toc(tBatch);
        est_total = elapsed * (NT / done_pts);
        eta = max(0, est_total - elapsed);
        fprintf('[GPU Capacity] batch %d/%d | pts %d | batch %.2fs | elapsed %.2fs | ETA ~%.2fs\n', ...
            b, total_batches, K, batch_t, elapsed, eta);
    end

else
    % -------- CPU parfor path WITH PROGRESS --------
    tx_flat_d = tx_flat(:);      % double for accuracy
    w_tx_d    = w_tx;            % double ok here
    lambda_d  = lambda;
    rx_off_d  = rx_offsets;      % double
    rho_d     = double(rho);
    Ns_d      = double(Ns);

    % linearize grid points (X fastest)
    [X0grid_d, Y0grid_d] = meshgrid(xv, yv);   % [Ny x Nx]
    X0v_d = reshape(X0grid_d.', [], 1);        % NT x 1
    Y0v_d = reshape(Y0grid_d.', [], 1);        % NT x 1

    C_lin = zeros(NT,1);
    p = gcp('nocreate'); if isempty(p), parpool; end

    % DataQueue progress (like your SVD section)
    qCap = parallel.pool.DataQueue; 
    afterEach(qCap, @capacityProgressListener);
    send(qCap, struct('init',true,'total',NT,'tStart',tCapStart));

    parfor t = 1:NT
        iterStart = tic;

        x0 = X0v_d(t); y0 = Y0v_d(t);

        % RX coordinates for this point
        xr = x0 + rx_off_d;     % 1 x N_r
        yr = y0;                % scalar

        % Element channel [num_tx x N_r]
        XE = reshape(tx_flat_d, [], 1);
        XR = reshape(xr, 1, []);
        R  = sqrt( (XR - XE).^2 + (yr - 0).^2 );
        H_elem = (1 ./ sqrt(4*pi*R)) .* exp(-1j * 2*pi .* R ./ lambda_d);

        % Effective subarray channel: H_eff [N_r x N_pa_tx]
        H_eff = (w_tx_d * H_elem).';   % (N_pa_tx x N_r)^T

        % Capacity via SVD (stable & fast)
        svals = svd(H_eff, 'econ');
        C_lin(t) = sum(log2(1 + (rho_d/Ns_d) * (svals.^2)));

        send(qCap, toc(iterStart));
    end

    % reshape back to [Ny x Nx] (X fastest)
    C_map = reshape(C_lin, [Nx, Ny]).';
end

fprintf('Capacity computation completed in %.2f s\n', toc(tCapStart));
%% ---------- Plot Capacity ----------
figure('Color','w');
imagesc(xv, yv, C_map);
if log_scale_axis, set(gca,'YDir','normal','YScale','log'); else, set(gca,'YDir','normal'); end
colormap(jet); cb=colorbar; cb.Label.String='b/s/Hz';
clim([min_capacity max_capacity]);
xlabel('Horizontal Offset (m)'); ylabel('Vertical RX–TX Separation (m)');
title(sprintf('SU-MIMO Capacity (N_r=%d, N_s=%d)', N_r, N_pa_tx));
hold on;
plot(poi_x, poi_y, 'ro', 'MarkerSize', 10, 'LineWidth', 2, ...
    'DisplayName','Focal Point');
legend('Location','northeastoutside');
y_overlay = min(vertical_offset_vals);
hold on;
plot(tx_flat(:), y_overlay*ones(numel(tx_flat),1), 'kp', ...
    'MarkerFaceColor','w','MarkerSize',8,'DisplayName','TX elements');

% Optional: show a quick table
fprintf('\n=== Per-Array SVD at Focal Point (x=%.2f m, y=%.2f m) ===\n', poi_x, poi_y);
for m = 1:N_pa_tx
    fprintf('Array %d: |E| = %.4g  (%.2f dB),  phase = %+6.1f°\n', ...
        m, sv_per_array(m), sv_per_array_dB(m), rad2deg(angle(E_focus(m))));
end
fprintf('Combined (||E||2): %.4g  (%.2f dB)\n\n', sv_joint, sv_joint_dB);


%% ---------- Helpers ----------
function progressListener(data)
    persistent completed total tStart
    if isstruct(data) && isfield(data,'init') && data.init
        completed = 0; total = data.total; tStart = data.tStart; return;
    end
    completed = completed + 1;
    fprintf('[CPU SVD] %.1f%%%% (%d/%d) | iter %.3fs | total %.1fs\n', ...
        100*completed/total, completed, total, data, toc(tStart));
end

function capacityProgressListener(data)
    persistent completed total tStart
    if isstruct(data) && isfield(data,'init') && data.init
        completed = 0; total = data.total; tStart = data.tStart; return;
    end
    completed = completed + 1;
    fprintf('[CPU Capacity] %.1f%%%% (%d/%d) | iter %.3fs | total %.1fs\n', ...
        100*completed/total, completed, total, data, toc(tStart));
end

function out = ternary(cond,a,b), if cond, out=a; else, out=b; end
end

function x = db2lin(dB), x = 10.^(dB/10); end

function C_bpsHz = su_mimo_capacity_at_point(x0, y0, tx_flat, w_tx, N_pa_tx, N_r, rx_spacing, lambda, B_Hz, Tx_power_W, T_kelvin, noiseFig_dB, k_B, useGPU,use_MRT_across_subarrays)
    % Build Rx antenna coords (ULA centered at x0)
    xr = x0 + ( (0:N_r-1) - (N_r-1)/2 ) * rx_spacing;
    yr = y0 * ones(1, N_r);

    % Per-element channels to each Rx antenna: H_elem_full [num_tx x N_r]
    num_tx = numel(tx_flat);
    XE = reshape(tx_flat(:), [num_tx,1]);
    YE = zeros(num_tx,1);

    XR = reshape(xr(:).', [1, N_r]);
    YR = reshape(yr(:).', [1, N_r]);

    DX = XR - XE;                            % [num_tx x N_r]
    DY = YR - YE;
    R  = sqrt(DX.^2 + DY.^2);
    H_elem_full = (1 ./ sqrt(4*pi*R)) .* exp(-1j * 2*pi .* R ./ lambda);   % [num_tx x N_r]

    % Combine within each subarray using your digital weights: W [N_pa_tx x num_tx]
    % Effective channel from subarrays to Rx antennas: H_eff [N_r x N_pa_tx]
    H_eff = (w_tx * H_elem_full).';   % (N_pa_tx x N_r)^T => [N_r x N_pa_tx]

    % Noise power with NF:
    N0B_W = k_B * T_kelvin * B_Hz * 10^(noiseFig_dB/10);
    rho   = Tx_power_W / N0B_W;       % total SNR (linear)

    if use_MRT_across_subarrays
        % 1-stream beamforming upper bound (not multiplexing): equal power across subarrays not used.
        g = vecnorm(H_eff,2,1);              % per-subarray gains if you were to combine? (not needed)
        % For MRT across subarrays (single stream), effective gain is largest singular value of H_eff:
        s = svd(H_eff);
        C_bpsHz = log2( 1 + rho * (s(1)^2) );   % all power to strongest mode
    else
        % Equal power across N_s = N_pa_tx independent streams:
        Ns = size(H_eff,2);
        C_bpsHz = real(log2(det( eye(N_r) + (rho/Ns) * (H_eff*H_eff') )));
    end

    
end