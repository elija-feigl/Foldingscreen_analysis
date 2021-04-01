%% ensure subfolders are in path
folder = fileparts(mfilename('fullpath'));
addpath(genpath(folder));

%% execute this to select a folder and analyze a folding screen
[pname] = uigetdir(pwd, 'Select an Directory with intial foling screen');
sigma_integrate_band = 1.0;

[name] = make_profiles(pname);
anotate_profiles(name, pname, sigma_integrate_band);
analyse_profiles(name, pname);
