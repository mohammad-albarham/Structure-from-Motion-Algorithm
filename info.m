function info(message, varargin)
% INFO - Custom function to display hierarchical information messages.
% Prints messages based on the global info_level and enableInfo settings.
%
% Usage:
%   info('Message')                        % Defaults to Level 1
%   info('Message', 1)                      % Explicit Level 1
%   info('Message with %d', 1, 10)          % Level 1 with variables
%   info('Message with %f', 2, 3.14159)     % Level 2 with variables
%
% Global Settings:
%   - enableInfo: Toggle info messages (true/false).
%   - info_level: Controls verbosity (1 = only Level 1, 2 = Level 1 and 2).

% Access global variables
global enableInfo;
global info_level;

% Set default level to 1 if not explicitly provided
if nargin < 2 || ~isnumeric(varargin{1}) % Check if level is missing
    level = 1;               % Default level
    args = varargin;         % Treat everything else as variables
else
    level = varargin{1};     % Extract level
    args = varargin(2:end);  % Remaining are variables
end

% Print messages based on levels
if enableInfo
    % Print if the current level is less than or equal to the global level
    if level <= info_level
        fprintf('[INFO][Level %d]: ', level);
        fprintf(message, args{:}); % Pass additional arguments
        fprintf('\n');
    end
end
end
