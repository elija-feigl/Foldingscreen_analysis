%% ensure subfolders are in path
folder = fileparts(mfilename('fullpath'));
addpath(genpath(folder));

%% analyse the profiles
[fname, pname] = uigetfile([pwd '/*.mat'], 'Select mat file');
name = fname(1:end-9);

fname = [name  '_data.mat'];
data = load([pname filesep fname]);

analyse_profiles(name, pname, data);