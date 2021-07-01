function data = compute_metrics(name, pname, data)
% step3: IFS analysis pipeline
% analyse the profiles
    fname = [name  '_data.mat'];
    
    % integrate aggregates, smear, and monomer bands
    data.gelInfo = integrate_species(data.gelInfo);
    data.gelInfo = metrics(data.gelInfo);
    
    % save
    disp(['Saving to ' pname filesep fname]);
    save([pname filesep fname], '-struct', 'data');
    disp('Done')
    
end