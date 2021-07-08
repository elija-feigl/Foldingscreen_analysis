%% ensure subfolders are in path
folder = fileparts(mfilename('fullpath'));
addpath(genpath(folder));

%% anotate the profiles by marking species
[fname, pname] = uigetfile([pwd '/*.mat'], 'Select mat file');
sigma_integrate_band = 1.0;
data = load([pname fname]); % load data for completion check
name = fname(1:end-9);


try
    if data.gelInfo.peaks_ok
        is_done = ~strcmp(questdlg('Pocket and monomers have already been selected. Redo?','Redo?' ,'No','Yes', 'Yes'),'Yes');
        if ~is_done
            data.gelInfo.peaks_ok = false; 
        end
    end
catch
    data.gelInfo.peaks_ok = false;      
end

if ~data.gelInfo.peaks_ok
    data = anotate_profiles(name, pname, data);
else
    disp('Nothing to do');
end
