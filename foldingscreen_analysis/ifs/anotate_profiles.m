function data = anotate_profiles(name, pname, data)
% step2: IFS analysis pipeline
% anotate the profiles by marking species
    fname = [name  '_data.mat'];
 
    try
        data.gelInfo.peaks_ok
    catch
        data.gelInfo.peaks_ok = false;      
    end

    % select aggregates and monomer bands
    [data.profileData, data.gelInfo] = select_species(data.profileData, data.gelData, data.gelInfo); %TODO: delete profileData


    % integrate aggregates, smear, and monomer bands
    [data.profileData, data.gelInfo] = integrate_species(data. profileData, data.gelInfo); %TODO: delete profileData
    
    % save
    disp(['Saving to ' pname filesep fname]);
    save([pname filesep fname], '-struct', 'data');
    disp('Done')

end
