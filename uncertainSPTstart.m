function uncertainSPTstart()
% add uncertainSPT to the matlab path

mfilepath=fileparts(mfilename('fullpath'));
addpath(genpath(mfilepath))
rmpath(genpath(fullfile(mfilepath,'.git')))

