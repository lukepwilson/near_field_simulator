function [first_eig_db, second_eig_db, eig_ratio_db, throughput, SNR_dB, tx_loc] = pa_nearfield_flexible_core(fc,dtrx,N_pa_rx,N_pa_tx,spread_rx,spread_tx,offset,beam_angle,null_angle)


% This function is for simulating multiple Phased Arrays placed in
% customizable configurations

% fc: center frequency, [constant]
% dtrx: distance between tx and rx, [constant]
% N_pa_rx: number of rx phased arrays, [constant]
% N_pa_tx: number of tx phased arrays, [constant]
% spread_rx: spread between the rx phased arrays, [constant]
% spread_tx: spread between the tx phased arrays, [constant]
% offset: horizontal offset between tx and rx arrays, [constant]
% beam_angle: beam angle for all phased arrays, [list] * WORK IN PROGRESS
% null_angle: null angle for all phased arrays, [list] * WORK IN PROGRESS


lambda = 3e8/fc;
d_int = lambda/2;

% ————————————————
% Create raw RX & TX placements (no centering yet)
% ————————————————
rx_phased_arrays_loc = zeros(N_pa_rx,8);
j = 0;
for i = 1:N_pa_rx
    rx_phased_arrays_loc(i,:) = (j:(j+7)) * d_int;
    j = j + 8;
end

tx_phased_arrays_loc = zeros(N_pa_tx,8);
j = 0;
for i = 1:N_pa_tx
    tx_phased_arrays_loc(i,:) = (j:(j+7)) * d_int;
    j = j + 8;
end


% ————————————————
% Apply spread and offset
% ————————————————
% RX: add spread_rx per subarray row, then the rx-offset from tx
rx_spread = rx_phased_arrays_loc + (0:N_pa_rx-1)' * spread_rx; 

% TX: add spread_tx per subarray row (no extra offset)
tx_spread = tx_phased_arrays_loc + (0:N_pa_tx-1)' * spread_tx;


% ————————————————
% NOW center *these* final positions about zero
% ————————————————
% RX
rx_flat = reshape(rx_spread,1,[]);                % all elements in one vector
rx_flat = rx_flat - mean(rx_flat);     % shift so avg = 0
rx_phased_arrays_loc_spread_offset = reshape(rx_flat, size(rx_spread)) + offset;

% TX
tx_flat = reshape(tx_spread,1,[]);
tx_flat = tx_flat - mean(tx_flat);
tx_phased_arrays_loc_spread = reshape(tx_flat, size(tx_spread));



% convert to a 2d coordinate matrix
rx_pa_loc2d = zeros([2 numel(rx_phased_arrays_loc_spread_offset)]);
tx_pa_loc2d = zeros([2 numel(tx_phased_arrays_loc_spread)]);


k = 1;
for j = 1:N_pa_rx
    for i=1:8
        rx_pa_loc2d(1,k) = rx_phased_arrays_loc_spread_offset(j,i);
        k=k+1;
    end
end


k = 1;
for j = 1:N_pa_tx
    for i=1:8
        tx_pa_loc2d(1,k) = tx_phased_arrays_loc_spread(j,i);
        k=k+1;
    end
end


rx_pa_loc2d(2,:) = dtrx*ones([1 numel(rx_pa_loc2d(2,:))]);
tx_pa_loc2d(2,:) = zeros([1 numel(tx_pa_loc2d(2,:))]);


% Build Channel Matrix
H = zeros(length(tx_pa_loc2d),length(rx_pa_loc2d));

for ii = 1:length(tx_pa_loc2d)
    for jj = 1:length(rx_pa_loc2d)
    
    x1 = tx_pa_loc2d(1,ii);
    x2 = rx_pa_loc2d(1,jj);
    y1 = tx_pa_loc2d(2,ii);
    y2 = rx_pa_loc2d(2,jj);

    d_ij = sqrt((y2-y1)^2+(x2-x1)^2);
    H(ii,jj) = exp(1j*2*pi*d_ij/lambda)*lambda/d_ij;
    
    end
end


eig_val = svd(H);
first_eig_db = db(eig_val(1));
second_eig_db = db(eig_val(2));


eig_ratio_db = first_eig_db - second_eig_db;

eig_ratio_db(eig_ratio_db>10) = 10;


eig_db = db(eig_val).';

Noise_power = -30;

SNR_dB = eig_db - Noise_power;



tx_power = 10;
SNR = db2pow(SNR_dB)';
alloc_power = waterfill(tx_power,1./SNR);
throughput = sum(log2(1+alloc_power'.*SNR'));

tx_loc = tx_pa_loc2d(1,:);



end

