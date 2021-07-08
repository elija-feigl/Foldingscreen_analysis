function analyse_profiles(name, pname, data)
% step3: IFS analysis pipeline
% analyse the profiles
    fname = [name  '_data.mat'];

    [data_tmp, cur_fig] = plot_analysis(data.gelInfo, data.gelData); 
    
    print(cur_fig, '-dpdf', [pname filesep name '_analysis.pdf']); %save figure
    data.foldingAnalysis = data_tmp;

    disp(['Saving to ' pname filesep fname])
    save([pname filesep fname], '-struct', 'data')
    disp('Done')
    
end

