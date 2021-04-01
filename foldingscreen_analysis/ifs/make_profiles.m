function [name] = make_profiles(pname)
% step1: IFS analysis pipeline
% process info of selected fodler
% process image 

    txt_files = dir([pname filesep '*.txt']);
    tif_files = dir([pname filesep '*.tif']);

    if length(tif_files)~=1 || length(txt_files)~=1
        disp(['Error. Found ' num2str(length(tif_files)) ' tif files and ' num2str(length(txt_files)) ' txt files in ' txt_files])
    else
        % get data
        txt = [txt_files(1).folder filesep txt_files(1).name];
        tif = [tif_files(1).folder filesep tif_files(1).name];
        [~, folder, ext] = fileparts(pname);
        name = strcat(folder, ext);
        name = split(name,'__');
        name = name{1};
        fname = [name  '_data.mat'];

        
        [gelData, gelInfo, profileData] = compute_profiles(pname, name, txt, tif);
        % save data
        disp('Saving data... please wait')
        save([pname filesep fname], 'gelData', 'gelInfo', 'profileData' )
        disp(['Data written to: ' [pname filesep fname]])
    end
end