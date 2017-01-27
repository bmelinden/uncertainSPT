function uncertainSPTstart()
% add uncertainSPT to the matlab path

dir0=fileparts(mfilename('fullpath'));% path to this file, even if called from other folder
addpath(dir0)
addpath(fullfile(dir0,'HMMcore'))
addpath(fullfile(dir0,'tools'))
disp('Added local uncertainSPT paths from')
disp(dir0)
disp('-----------------------')

