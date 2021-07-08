function data = anotate_profiles(name, pname, data)
% step2: IFS analysis pipeline
% anotate the profiles by marking species
    fname = [name  '_data.mat'];
 
    try
        data.gelInfo.peaks_ok
    catch
        data.gelInfo.peaks_ok = false;      
    end
    
    data.gelInfo.sigma_integrate = 1.0;
    
    while ~data.gelInfo.peaks_ok
        % select aggregates and monomer bands
        data.gelInfo = select_species(data.gelData, data.gelInfo);
    end
    
    % save
    disp(['Saving to ' pname filesep fname]);
    save([pname filesep fname], '-struct', 'data');
    disp('Done')

end
