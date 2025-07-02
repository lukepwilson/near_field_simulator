% PHASED ARRAY FLEXIBLE SIMULATION WRAPPER (4-panel)
% Run the script and interact with the inputs to update the four heatmaps
% Author: Luke Wilson
% Date: 07/2025

clear; clc; close all;

%% INITIAL PARAMETERS
P = struct();
P.fc                     = 29e9;    % center frequency (Hz)
P.N_pa_rx                = 2;       % number of RX arrays
P.N_pa_tx                = 2;       % number of TX arrays
P.spread_rx              = 0.10;    % RX spread (m)
P.spread_tx              = 0.05;    % TX spread (m)

P.vertical_offset_vals   = linspace(0.001, 10, 150);   % Vertical offset (m)
P.horizontal_offset_vals = linspace(-2, 2, 150);       % Horizontal offset (m)

%% FIGURE & TILES SETUP
fig = figure('Name','Flexible Core: 4 Metrics','NumberTitle','off', ...
             'Position',[200 200 1000 800]);

t = tiledlayout(fig,2,2, ...
      'TileSpacing','compact', ...
      'Padding','compact', ...
      'Position',[0.05 0.30 0.90 0.65]);

ax1 = nexttile(t); title(ax1,'Capacity (bits/s/Hz)','FontSize',14);
ax2 = nexttile(t); title(ax2,'Eig Ratio (dB)','FontSize',14);
ax3 = nexttile(t); title(ax3,'Eig 1 (dB)','FontSize',14);
ax4 = nexttile(t); title(ax4,'Eig 2 (dB)','FontSize',14);

infoBox = annotation(fig,'textbox',[0.70 0.30 0.25 0.08], ...
    'String',{}, 'EdgeColor','none', 'HorizontalAlignment','right', 'FontSize',12);

setappdata(fig,'P',P);
setappdata(fig,'Axes',[ax1, ax2, ax3, ax4]);
setappdata(fig,'InfoBox',infoBox);

%% UI CONTROLS (bottom 30%)
uiFontSize = 12;
% # TX arrays
uicontrol(fig,'Style','edit','Units','normalized',...    
    'Position',[0.05 0.25 0.10 0.04],...                   
    'String',num2str(P.N_pa_tx),...                       
    'FontSize',uiFontSize,...
    'Callback',@(e,~) editCallback(fig,'N_pa_tx',e));
uicontrol(fig,'Style','text','Units','normalized',...    
    'Position',[0.17 0.25 0.15 0.04],...                   
    'String','# TX arrays','FontSize',uiFontSize,...
    'HorizontalAlignment','left');

% # RX arrays
uicontrol(fig,'Style','edit','Units','normalized',...    
    'Position',[0.05 0.19 0.10 0.04],...                   
    'String',num2str(P.N_pa_rx),...                       
    'FontSize',uiFontSize,...
    'Callback',@(e,~) editCallback(fig,'N_pa_rx',e));
uicontrol(fig,'Style','text','Units','normalized',...    
    'Position',[0.17 0.19 0.15 0.04],...                   
    'String','# RX arrays','FontSize',uiFontSize,...
    'HorizontalAlignment','left');

% Center Frequency slider
uicontrol(fig,'Style','slider','Units','normalized',...   
    'Position',[0.05 0.14 0.80 0.04],...                   
    'Min',25e9,'Max',30e9,'Value',P.fc,...
    'Callback',@(s,~) sliderCallback(fig,'fc',get(s,'Value')));
uicontrol(fig,'Style','text','Units','normalized',...    
    'Position',[0.87 0.14 0.10 0.04],...                   
    'String','Center Frequency (Hz)','FontSize',uiFontSize,...
    'HorizontalAlignment','left');

% RX Spread slider
uicontrol(fig,'Style','slider','Units','normalized',...   
    'Position',[0.05 0.09 0.80 0.04],...                   
    'Min',0,'Max',1,'Value',P.spread_rx,...
    'Callback',@(s,~) sliderCallback(fig,'spread_rx',get(s,'Value')));
uicontrol(fig,'Style','text','Units','normalized',...    
    'Position',[0.87 0.09 0.10 0.04],...                   
    'String','Spread RX (m)','FontSize',uiFontSize,...
    'HorizontalAlignment','left');

% TX Spread slider
uicontrol(fig,'Style','slider','Units','normalized',...   
    'Position',[0.05 0.04 0.80 0.04],...                   
    'Min',0,'Max',1,'Value',P.spread_tx,...
    'Callback',@(s,~) sliderCallback(fig,'spread_tx',get(s,'Value')));
uicontrol(fig,'Style','text','Units','normalized',...    
    'Position',[0.87 0.04 0.10 0.04],...                   
    'String','Spread TX (m)','FontSize',uiFontSize,...
    'HorizontalAlignment','left');

%% INITIAL DRAW
drawAll(fig);

%% CALLBACK: Slider change
function sliderCallback(fig,field,val)
    P = getappdata(fig,'P');
    P.(field) = val;
    setappdata(fig,'P',P);
    drawAll(fig);
end

%% CALLBACK: Editbox change
function editCallback(fig,field,src)
    P = getappdata(fig,'P');
    v = round(str2double(get(src,'String')));
    if isnan(v) || v < 1, v = P.(field); end
    set(src,'String',num2str(v));
    P.(field) = v;
    setappdata(fig,'P',P);
    drawAll(fig);
end

%% DRAW ALL 4 PANELS
function drawAll(fig)
    P   = getappdata(fig,'P');
    axs = getappdata(fig,'Axes');
    infoBox = getappdata(fig,'InfoBox');
    M = numel(P.vertical_offset_vals);
    N = numel(P.horizontal_offset_vals);

    cap_map    = nan(M,N);
    ratio_map  = nan(M,N);
    eig1_map   = nan(M,N);
    eig2_map   = nan(M,N);

    for i = 1:M
      for j = 1:N
        [first_eig_db, second_eig_db, eig_ratio_db, throughput, ~, tx_locs] = ...
          pa_nearfield_flexible_core( ...
            P.fc, ...
            P.vertical_offset_vals(i), ...    % vertical offset
            P.N_pa_rx, P.N_pa_tx, ...
            P.spread_rx, P.spread_tx, ...
            P.horizontal_offset_vals(j), ...  % horizontal offset
            [], [] );

        eig1_map(i,j)  = first_eig_db;
        eig2_map(i,j)  = second_eig_db;
        ratio_map(i,j) = eig_ratio_db;

        % throughput already = spectral efficiency [bits/s/Hz]
        cap_map(i,j)   = throughput;
      end
    end

    all_eigs = [eig1_map(:); eig2_map(:)];
    eig_min = min(all_eigs);
    eig_max = max(all_eigs);

    maps   = {cap_map, ratio_map, eig1_map, eig2_map};
    labels = {'Capacity','Eig Ratio (dB)','Eig 1 (dB)','Eig 2 (dB)'};

    for k = 1:4
      ax = axs(k);
      cla(ax);
      imagesc(ax, P.horizontal_offset_vals, P.vertical_offset_vals, maps{k});
      set(ax,'YDir','normal','FontSize',12);
      xlabel(ax,'Horizontal Offset (m)','FontSize',12);
      ylabel(ax,'Vertical Offset (m)','FontSize',12);
      title(ax,labels{k},'FontSize',14);

      if k==3 || k==4
        clim(ax,[eig_min eig_max]);
      end

      cb = colorbar(ax);
      if k==1
        cb.Label.String = 'Capacity (bits/s/Hz)';
      else
        cb.Label.String = labels{k};
      end
      cb.Label.FontSize = 12;

      hold(ax,'on');
      plot(ax, tx_locs, zeros(size(tx_locs)), 'kp','MarkerFaceColor','w');
      hold(ax,'off');
    end

    infoStr = sprintf('fc: %.2e Hz\nSpread RX: %.2f m\nSpread TX: %.2f m', ...
        P.fc, P.spread_rx, P.spread_tx);
    set(infoBox,'String',infoStr);
end
