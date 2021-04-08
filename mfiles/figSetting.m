function figSetting(varargin)
%
% Modification: now don't change the root definitions
% Modification 27/10/2020 : added several pre-sets and the text size

switch nargin
    case 0
        width    = 12;
        height   = 10;
        textSize = 12;
    case 1
        picStyle = varargin{1};
        switch picStyle
            case 'paper'   %Single column for IEEE papers
                %External size in cm (double column IEEE paper)
                width = 8.58;
                height = 4;
                % Text size. 8pt as figure caption size.
                textSize = 8;
            case 'paperDouble' %Single column, larger height (e.g. double plots
                %External size in cm (double column IEEE paper)
                width = 19;
                height = 6;
                % Text size. 8pt as figure caption size.
                textSize = 8;
            case 'large'  %Format Single Column for thesis
                width = 15;
                height = 7.5;
                textSize = 10;
            case 'a4vert'   %vertical A4 paper
                width = 21;
                height = 29.7;
                textSize = 12;
            case 'a4horiz'   %horizontal A4 paper
                width = 29.7;
                height = 21;
                textSize = 12;
            otherwise
                warning('off','backtrace')
                warning('Pre-set not known. Default settings applied.')
                warning('on','backtrace')
                width    = 12;
                height   = 10;
                textSize = 12;
        end
    case 2
        width    = varargin{1};
        height   = varargin{2};
        textSize = 12;
    case 3
        width    = varargin{1};
        height   = varargin{2};
        textSize = varargin{3};
    otherwise
        warning('off','backtrace')
        warning('Error in the number of inputs: default settings imposed')
        warning('on','backtrace')
        width    = 12;
        height   = 10;
        textSize = 12;
end

% Definitions for plotting figures

set(gcf,'defaultTextInterpreter','Latex');
set(gcf,'defaultLegendInterpreter','Latex');
set(gcf,'defaultAxesTickLabelInterpreter','Latex');

set(gcf,'defaultAxesLineWidth',1);
set(gcf,'defaultLineLineWidth',1.5);

set(gcf,'defaultAxesGridLineStyle',':');
set(gcf,'defaultAxesXGrid','on');
set(gcf,'defaultAxesYGrid','on');
set(gcf,'defaultAxesZGrid','on');
set(gcf,'defaultAxesXColor',0*[1 1 1]);
set(gcf,'defaultAxesYColor',0*[1 1 1]);
set(gcf,'defaultAxesZColor',0*[1 1 1]);
%set(gcf,'defaultAxesGridAlpha',1);
%set(gcf,'defaultAxesLayer','top');
set(gcf,'defaultAxesBox','on');
set(gcf,'defaultAxesNextPlot','add');

set(gcf,'defaultAxesFontSize',textSize);
set(gcf,'defaultTextFontSize',textSize);
set(gcf,'defaultLegendFontSize',textSize);
set(gcf,'defaultAxesFontSizeMode','manual');

set(gcf,'defaultTextFontSizeMode','manual');
set(gcf,'defaultLegendFontSizeMode','manual');
set(gcf,'defaultAxesLabelFontSizeMultiplier',1);
set(gcf,'defaultAxesTitleFontSizeMultiplier',1);

set(gcf,'defaultAxesFontName','Times');
set(gcf,'defaultTextFontName','Times');

screenPos=get(groot,'ScreenSize')/get(groot,'ScreenPixelsPerInch')*2.54; % cm
figPos(1)=screenPos(3)/2-width/2;
figPos(2)=screenPos(4)/2-height/2;
figPos(3)=width;
figPos(4)=height;

colors{1} = [0.0 0.0 1.0];
colors{2} = [1.0 0.0 0.0];
colors{3} = [0.0 0.8 0.0];
colors{4} = [1.0 0.5 0.0];
colors{5} = [0.0 0.8 0.8];
colors{6} = [0.8 0.0 0.8];
colors{7} = [1.0 0.8 0.0];

set(gcf,'defaultAxesColorOrder',[colors{1};colors{2};colors{3};colors{4};colors{5};colors{6};colors{7}]);

% set(gcf,'Renderer','painters');
set(gcf,'Units','centimeters');
set(gcf,'Position',figPos);
set(gcf,'Color',[1 1 1]);
set(gcf,...
    'PaperUnits','centimeter',...
    'PaperPosition',[0 0 width height],...
    'PaperType','<custom>',...
    'PaperSize',[width height])
