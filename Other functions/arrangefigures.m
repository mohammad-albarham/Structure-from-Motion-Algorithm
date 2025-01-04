function arrangefigures(varargin)
% arrangeFigures automatically aranges open figures in a grid on the
% designated screen monitor. Options are:
%
% arrangefigures
% arrangefigures(monitor)
% arrangefigures(ontop)
% arrangefigures(monitor, fRatio)
% arrangefigures(monitor, [w h])
% arrangefigures(monitor, fRatio, animate)
% 
%
% Possible settings are:
%
% monitor:
% - default: 1
% - 0 means all monitors
% - scalar input means specified monitor
% - array input means multiple chosen monitors. Monitors must be successive
% ontop:
% - if there is only one input argument and the input argument is of type
%   logical and is set to true, then all active figures are placed on top.
%   Everything else is left unchanged.
% fRatio:
% - default: 3/4
% - set the figure ratio in height/width (default: 3/4)
% [h w]:
% - default: -> default is fRatio
% - set a fixed width (w) and height (h) in pixels
% animate:
% - default: false
% - boolean flag, added for fun, don't expect too much... Can also be of
%   practical use when all figures needs to be place on top.
%
% (c) M. van Dijk

if nargin > 3
    warning('Incorrect number of input arguments. Leaving figures unchanged.')
    return
elseif nargin == 0
    monitor = 1;
    fRatio = 3/4;
    animate = false;
elseif nargin == 1
    monitor = varargin{1};
    fRatio = 3/4;
    animate = false;
elseif nargin == 2
    monitor = varargin{1};
    fRatio = varargin{2};
    animate = false;
elseif nargin == 3
    monitor = varargin{1};
    fRatio = varargin{2};
    animate = varargin{3};
end

if numel(fRatio) == 2
    fh = fRatio(2);
    fw = fRatio(1);
    fRatio = fh/fw;
end

% monitor = 2;
% fRatio = 0.25;
% animate = false;
spacer = 30;

% get figure handles
figureHandles = sort(get(0, 'Children'));
NFigures = numel(figureHandles);

% if first argument is logical, only place all figures on top
if islogical(monitor)
    for figureHandle = figureHandles'
        figure(figureHandle);
    end
    return
end

% get number of monitors and the resolutions
MonitorExtends = get(0, 'MonitorPositions');
NMonitors = size(MonitorExtends, 1);
if isscalar(monitor)
    if monitor == 0
        monitor = 1:NMonitors;
    end
end

% check if monitor-array exceeds number of monitors
if min(monitor) > NMonitors
    error(['Can''t find assigned monitor: NMonitors = ' num2str(NMonitors) ', Set display monitor = [' num2str((monitor)) ']'])
end
if min(monitor) < 1
    warning(['Can''t find assigned monitor: NMonitors = ' num2str(NMonitors) ', Set display monitor = [' num2str((monitor)) ']. Skipping monitors < 0'])
    monitor = monitor(monitor > 0);    
end
if max(monitor) > NMonitors
    warning(['Can''t find assigned monitor: NMonitors = ' num2str(NMonitors) ', Set display monitor = [' num2str((monitor)) ']. Skipping monitors > ' num2str(NMonitors)])
    monitor = monitor(monitor <= NMonitors);
end
if isempty(monitor)
    error('Assign a correct monitor')
end

UsedMonitorExtends = MonitorExtends(monitor, :);

taskbarOffset = 100;
monitorArea = [min(UsedMonitorExtends(:, 1:2), [], 1) max(UsedMonitorExtends(:, 3:4), [], 1)];
x = monitorArea(1);
y = monitorArea(2);
w = monitorArea(3) - x;
h = monitorArea(4) - taskbarOffset - y; % 100 offset for taskbar

% determine size figures and grid
hSpacer = 75;
if ~exist('fw', 'var')
    figureArea = 0;
    for NRows = 1:NFigures
        NCols = ceil(NFigures / NRows);
        fh = (h - 2*spacer) / NRows - spacer - hSpacer;
        fw = fh / fRatio;
        if fw > ((w - 2*spacer) / NCols - spacer)
            fw = (w - 2*spacer) / NCols - spacer;
            fh = fw * fRatio;
        end

        if figureArea < fh*fw
            setNRows = NRows;
            setNCols = NCols;
            figureArea = fh*fw;
        end
    end
    NRows = setNRows;
    NCols = setNCols;
    fh = (h - 2*spacer) / NRows - spacer - hSpacer;
    fw = (w - 2*spacer) / NCols - spacer;
    if (fh/fw > fRatio)
        fh = fw * fRatio;
    else
        fw = fh / fRatio;
    end
else
    NCols = floor((w - 2*spacer) / (fw - spacer));
    NRows = floor((h - 2*spacer) / (fh - spacer));
    if (NCols * NRows < NFigures)
        warning(['Cannot display ' num2str(NFigures) ' figures using [w h] = [' num2str(w) ' ' num2str(h) '] on monitor ' num2str(monitor) '. Leaving figures unchanged.']);
        return
    end
end

fw = floor(repmat(fw, [1 NFigures]));
fh = floor(repmat(fh, [1 NFigures]));
fx = repmat(spacer+(0:NCols-1)*(fw(1)+spacer), [1, NRows]) + x;
fy = (h - fh(1)) - sort(reshape(repmat(spacer+(0:NRows-1)*(fh(1)+spacer+hSpacer), [NCols, 1]), 1, NCols*NRows), 'ascend') + y;


if animate
    nsteps = 10;
    posarr = zeros(NFigures, 4, nsteps);
    for n = 1:NFigures
        posarr(n, :, 1) = get(figureHandles(n), 'position');
        posarr(n, 1, :) = linspace(posarr(n, 1, 1), fx(n), nsteps);
        posarr(n, 2, :) = linspace(posarr(n, 2, 1), fy(n), nsteps);
        posarr(n, 3, :) = linspace(posarr(n, 3, 1), fw(n), nsteps);
        posarr(n, 4, :) = linspace(posarr(n, 4, 1), fh(n), nsteps);
        figure(figureHandles(n));
    end
    for step = 1:nsteps
        for n = 1:NFigures
            set(figureHandles(n), 'position', squeeze(posarr(n, :, step)));
        end
        pause(0.007);
    end
else
    for n = 1:NFigures
        set(figureHandles(n), 'position', [fx(n) fy(n) fw(n) fh(n)]);
    end
end



