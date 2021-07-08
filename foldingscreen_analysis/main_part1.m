 %% ensure subfolders are in path
folder = fileparts(mfilename('fullpath'));
addpath(genpath(folder));

%% execute this to select a folder and analyze a folding screen
[pname] = uigetdir(pwd, 'Select an Directory with intial folding screen');
make_profiles(pname)