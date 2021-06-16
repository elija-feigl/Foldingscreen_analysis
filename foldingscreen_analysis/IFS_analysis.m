%% ensure subfolders are in path
folder = fileparts(mfilename('fullpath'));
addpath(genpath(folder));

%% execute this to select a folder and analyze a folding screen
[pname] = uigetdir(pwd, 'Select an Directory with intial folding screen');
sigma_integrate_band = 1.0;

% Step 1
[name, data] = make_profiles(pname);

% Step 2
data = anotate_profiles(name, pname, data);

% Step 3
analyse_profiles(name, pname, data);
