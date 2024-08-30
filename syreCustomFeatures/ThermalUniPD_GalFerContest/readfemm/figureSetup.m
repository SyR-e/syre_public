function figureSetup(h0,varargin)

% FIGURESETUP configure matlab figuers with standard template
%
% USE:
% figureSetup(h)
%
% INPUT:
% 'h': figure handle
% 'varargin': variable inputs
%   * 'FontSize': font size (default: 14 pts)
%   * 'LineWidth': line width (default: 2 pts)
%   * 'MarkerSize': marker size (default: 8 pts)
%   * 'Interpreter': interpreter. Available: 'tex', 'latex' (default: 'tex')
%
% OUTPUT:
%
% VERSION:
% Date: 12.04.2012
% Copyright(C) 2012-2023: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 25.05.2012: added white background
% 12.09.2012: added variable inputs
% 21.07.2014: added ZOOM as optional input
% 20.07.2016: refurbished routine working also with subplot
% 12.01.2011: added BOX as optional input
% 15.05.2023: removed GRID and BOX fro optional inputs
% 08.06.2023: added INTERPRETER as optional input
% 31.07.2023: fixed bug with SUBPLOTS

% set defaults
FontSize = 16;
LineWidth = 2;
MarkerSize = 8;
Zoom = 1;
Interpreter = 'tex';

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'FontSize')
        FontSize = varargin{i+1};
    elseif strcmpi(vin,'LineWidth')
        LineWidth = varargin{i+1};
    elseif strcmpi(vin,'MarkerSize')
        MarkerSize = varargin{i+1};
    elseif strcmpi(vin,'Zoom')
        Zoom = varargin{i+1};
    elseif strcmpi(vin,'Interpreter')
        Interpreter = varargin{i+1};
    end
end

% get figure handle
figure(h0);

% axis equal;
if length(axis(gca)) == 4
    grid = 'on';
    box = 'on';
else
    grid = 'off';
    box = 'off';
end

% zoom
zoom(Zoom);

% get children
h1 = get(h0,'children');

% loop over children
for i = 1:length(h1)
    if ismethod(h1(i),'Axes')
        % fontsize
        h1(i).FontSize = FontSize;
        % box on
        h1(i).Box = box;
        % grid off
        h1(i).XGrid = grid;
        h1(i).YGrid = grid;
        h1(i).ZGrid = grid;
        % fix linewidth
        set(findobj(h1(i),'Type','line'),'LineWidth',LineWidth);
        % marker size
        set(findobj(h1(i),'Type','line'),'Markersize',MarkerSize)
        % fix free text
        set(findobj(h1(i),'Type','text'),'FontSize',FontSize)
    elseif ismethod(h1(i),'Legend')
        h1(i).FontSize = .9*FontSize;
    end
end

% set interpreter
set(gca,'TickLabelInterpreter',Interpreter);
set(findobj(gcf, 'Type', 'Legend'),'Interpreter',Interpreter);

hxlab = get(findobj(gcf,'type','axe'),'xlabel');
if isscalar(hxlab)
    set(hxlab,'Interpreter',Interpreter);
else
    for i = 1:length(hxlab)
        set(hxlab{i},'Interpreter',Interpreter);
    end
end
hylab = get(findobj(gcf,'type','axe'),'ylabel');
if isscalar(hylab)
    set(hylab,'Interpreter',Interpreter);
else
    for i = 1:length(hylab)
        set(hylab{i},'Interpreter',Interpreter);
    end
end

% set(get(findobj(gcf,'type','axe'),'ylabel'),'Interpreter',Interpreter);

% set white background
set(h0,'Color',[1 1 1]);

