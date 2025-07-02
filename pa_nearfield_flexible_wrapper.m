% PHASED ARRAY FLEXIBLE SIMULATION WRAPPER (4-panel)
% Author: Luke Wilson | UI Optimized: 07/2025

clear; clc; close all;

%% INITIAL PARAMETERS
P = struct();
P.fc               = 29e9;
P.N_pa_rx          = 1;
P.N_pa_tx          = 1;
P.spread_rx        = 0.05;
P.spread_tx        = 0.05;
P.horizontal_range = 0.1;
P.vertical_range   = 0.4;
P.vertical_offset_vals   = linspace(0.001, P.vertical_range, 150);
P.horizontal_offset_vals = linspace(-P.horizontal_range, P.horizontal_range, 150);

%% FIGURE & TILES SETUP
fig = figure('Name','Flexible Core: 4 Metrics','NumberTitle','off', ...
             'Position',[200 200 1000 900]);

t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact', ...
                'Position',[0.05 0.25 0.90 0.7]);  % Enlarged and moved up

ax1 = nexttile(t); title(ax1,'Capacity (bits/s/Hz)','FontSize',14);
ax2 = nexttile(t); title(ax2,'Eig Ratio (dB)','FontSize',14);
ax3 = nexttile(t); title(ax3,'Eig 1 (dB)','FontSize',14);
ax4 = nexttile(t); title(ax4,'Eig 2 (dB)','FontSize',14);

infoBox = annotation(fig,'textbox',[0.70 0.45 0.25 0.08], ...
    'String',{}, 'EdgeColor','none', 'HorizontalAlignment','right', 'FontSize',12);

setappdata(fig,'P',P);
setappdata(fig,'Axes',[ax1, ax2, ax3, ax4]);
setappdata(fig,'InfoBox',infoBox);

%% UI LAYOUT (2 columns, compact)
uiFontSize = 12;
halfW = 0.34; spacing = 0.10; labelW = 0.12;
left1 = 0.05; left2 = left1 + halfW + spacing;
fullW = 0.79;

% Updated row vertical positions (lowered)
y1 = 0.15;  % Row 1: TX / RX arrays
y2 = 0.10;  % Row 2: Center Frequency
y3 = 0.05;  % Row 3: Spread TX / RX
y4 = 0.02;  % Row 4: Horizontal / Vertical offset range sliders

% --- Row 1: TX / RX Arrays ---
uicontrol(fig,'Style','edit','Units','normalized',...
    'Position',[left1 y1 0.07 0.04],'String',num2str(P.N_pa_tx),'FontSize',uiFontSize,...
    'Callback',@(e,~) editCallback(fig,'N_pa_tx',e));
uicontrol(fig,'Style','text','Units','normalized',...
    'Position',[left1+0.08 y1 labelW 0.04],'String','# TX arrays','FontSize',uiFontSize,...
    'HorizontalAlignment','left');

uicontrol(fig,'Style','edit','Units','normalized',...
    'Position',[left2 y1 0.07 0.04],'String',num2str(P.N_pa_rx),'FontSize',uiFontSize,...
    'Callback',@(e,~) editCallback(fig,'N_pa_rx',e));
uicontrol(fig,'Style','text','Units','normalized',...
    'Position',[left2+0.08 y1 labelW 0.04],'String','# RX arrays','FontSize',uiFontSize,...
    'HorizontalAlignment','left');

% --- Row 2: Center Frequency (full width) ---
uicontrol(fig,'Style','slider','Units','normalized',...
    'Position',[left1 y2 fullW 0.035],'Min',25e9,'Max',30e9,'Value',P.fc,...
    'Callback',@(s,~) sliderCallback(fig,'fc',get(s,'Value')));
uicontrol(fig,'Style','text','Units','normalized',...
    'Position',[left1+fullW+0.01 y2 labelW 0.035],'String','Center Frequency (Hz)', ...
    'FontSize',uiFontSize,'HorizontalAlignment','left');

% --- Row 3: Spread TX / Spread RX ---
uicontrol(fig,'Style','slider','Units','normalized',...
    'Position',[left1 y3 halfW 0.035],'Min',0,'Max',1,'Value',P.spread_tx,...
    'Callback',@(s,~) sliderCallback(fig,'spread_tx',get(s,'Value')));
uicontrol(fig,'Style','text','Units','normalized',...
    'Position',[left1+halfW+0.005 y3 labelW 0.035],'String','Spread TX (m)', ...
    'FontSize',uiFontSize,'HorizontalAlignment','left');

uicontrol(fig,'Style','slider','Units','normalized',...
    'Position',[left2 y3 halfW 0.035],'Min',0,'Max',1,'Value',P.spread_rx,...
    'Callback',@(s,~) sliderCallback(fig,'spread_rx',get(s,'Value')));
uicontrol(fig,'Style','text','Units','normalized',...
    'Position',[left2+halfW+0.005 y3 labelW 0.035],'String','Spread RX (m)', ...
    'FontSize',uiFontSize,'HorizontalAlignment','left');

% --- Row 4: Horizontal / Vertical Offset Ranges ---
uicontrol(fig,'Style','slider','Units','normalized',...
    'Position',[left1 y4 halfW 0.035],'Min',0.01,'Max',2,'Value',P.horizontal_range,...
    'Callback',@(s,~) rangeSliderCallback(fig,'horizontal_range',get(s,'Value')));
uicontrol(fig,'Style','text','Units','normalized',...
    'Position',[left1+halfW+0.005 y4 labelW 0.035],'String','Horz Offset (±m)', ...
    'FontSize',uiFontSize,'HorizontalAlignment','left');

uicontrol(fig,'Style','slider','Units','normalized',...
    'Position',[left2 y4 halfW 0.035],'Min',0.01,'Max',10,'Value',P.vertical_range,...
    'Callback',@(s,~) rangeSliderCallback(fig,'vertical_range',get(s,'Value')));
uicontrol(fig,'Style','text','Units','normalized',...
    'Position',[left2+halfW+0.005 y4 labelW 0.035],'String','Vert Offset Max (m)', ...
    'FontSize',uiFontSize,'HorizontalAlignment','left');

%% INITIAL DRAW
drawAll(fig);

%% CALLBACK FUNCTIONS
function sliderCallback(fig,field,val)
    P = getappdata(fig,'P');
    P.(field) = val;
    setappdata(fig,'P',P);
    drawAll(fig);
end

function editCallback(fig,field,src)
    P = getappdata(fig,'P');
    v = round(str2double(get(src,'String')));
    if isnan(v) || v < 1, v = P.(field); end
    set(src,'String',num2str(v));
    P.(field) = v;
    setappdata(fig,'P',P);
    drawAll(fig);
end

function rangeSliderCallback(fig, field, val)
    P = getappdata(fig,'P');
    P.(field) = val;
    if strcmp(field,'horizontal_range')
        P.horizontal_offset_vals = linspace(-val, val, 150);
    elseif strcmp(field,'vertical_range')
        P.vertical_offset_vals = linspace(0.001, val, 150);
    end
    setappdata(fig,'P',P);
    drawAll(fig);
end

%% DRAW FUNCTION
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
        [eig1, eig2, eigRatio, capacity, ~, tx_locs] = ...
          pa_nearfield_flexible_core( ...
            P.fc, P.vertical_offset_vals(i), ...
            P.N_pa_rx, P.N_pa_tx, ...
            P.spread_rx, P.spread_tx, ...
            P.horizontal_offset_vals(j), [], [] );

        eig1_map(i,j)  = eig1;
        eig2_map(i,j)  = eig2;
        ratio_map(i,j) = eigRatio;
        cap_map(i,j)   = capacity;
      end
    end

    all_eigs = [eig1_map(:); eig2_map(:)];
    eig_min = min(all_eigs);
    eig_max = max(all_eigs);
    maps   = {cap_map, ratio_map, eig1_map, eig2_map};
    labels = {'Capacity','Eig Ratio (dB)','Eig 1 (dB)','Eig 2 (dB)'};

    for k = 1:4
        ax = axs(k); cla(ax);
        imagesc(ax, P.horizontal_offset_vals, P.vertical_offset_vals, maps{k});
        set(ax,'YDir','normal','FontSize',12);
        xlabel(ax,'Horizontal Offset (m)','FontSize',12);
        ylabel(ax,'Vertical Offset (m)','FontSize',12);
        title(ax,labels{k},'FontSize',14);

        if k==3 || k==4
            clim(ax,[eig_min eig_max]);
        end
        cb = colorbar(ax);
        cb.Label.String = labels{k}; cb.Label.FontSize = 12;

        hold(ax,'on');
        plot(ax, tx_locs, zeros(size(tx_locs)), 'kp','MarkerFaceColor','w');
        hold(ax,'off');
    end

    infoStr = sprintf(['fc: %.2e Hz\nSpread RX: %.2f m\nSpread TX: %.2f m\n' ...
                       'Horz Range: ±%.2f m\nVert Range: %.2f m'], ...
        P.fc, P.spread_rx, P.spread_tx, P.horizontal_range, P.vertical_range);
    set(infoBox,'String',infoStr);
end