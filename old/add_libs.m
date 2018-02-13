currentDir = mfilename('fullpath');

% normalise the slash
currentDir = strrep(currentDir, '\', '/');

% parse the directory only
kkkkk = strfind(currentDir, '/');
currentDir = currentDir(1:kkkkk(end));

% add all paths
if ~strcmp(currentDir(end), '/') currentDir = strcat(currentDir, '/'); end
addpath(genpath(currentDir));
% addpath(currentDir);
% addpath(strcat(currentDir, 'lib'));
% addpath(strcat(currentDir, 'lib/minFunc_2012'));
% addpath(strcat(currentDir, 'lib/minFunc_2012/autoDif'));
