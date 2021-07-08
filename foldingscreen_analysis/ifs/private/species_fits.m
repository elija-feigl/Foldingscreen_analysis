        
%% Calculate fits
function gelInfo = species_fits(gelInfo, gelData)
    %% Get band fits
    function [fits, fit_range] = get_band_fit(species, pos)
        xpos = species.select_positions(pos) - fix(species.height/2);
        fits = zeros(1, 3);
        
        upper = xpos + species.height;
        val = species.fullprofiles{pos};
        
        if xpos < pocket_position
            fits(1) = 0.0;
            fits(2) = pocket_position;
            fits(3) = 0.0;
            fit_range = pocket_position;
        else
            y = val(xpos:upper);
            x = double(xpos:upper);
            options = fitoptions('gauss1', 'Lower', [0 xpos 0], 'Upper', [Inf upper Inf]);
            
            fit_range = [options.lower(2) options.upper(2)];
            fits = coeffvalues(fit(x', y, 'gauss1', options));
        end
    end

    %% Get staple fists
    function fits = get_staple_fit(species, pos, height)
        if species.type == "ladder"
            fits = zeros(1,3);
        elseif species.type == "scaffold"
            fits = zeros(1,3);
        else
            xpos = species.staple.select_positions(pos) - fix(species.height/2);
            if xpos < pocket_position
                fits = [0.0 pocket_position 0.0];
            else
                upper = xpos + height;
                val = species.fullprofiles{pos};   

                if upper > length(val)
                    upper = length(val);
                end

                y = val(xpos:upper);
                x = double(xpos:upper);
                options = fitoptions('gauss1', 'Lower', [0 xpos 0], 'Upper', [Inf upper Inf]);

                fits = coeffvalues(fit(x', y, 'gauss1', options));
            end
        end
    end

    % TODO: hardcoding these might not be ideal, as pixelsize can depend
    %       on the scanner used. find better solution!
    gelInfo.species.ladder.height = 50;
    gelInfo.species.scaffold.height = 50;
    gelInfo.species.mono.height = 50;
    stapleHeight = 200; 

    %% Pocket, Ladder, Scaffold, Monomers and staple fits
    
    % compute pocket sum profiles and fit it with gaussian.
    % position and width of pocket is always the same -> sum
    
    pocket_position = double(gelInfo.pocket.bounding_box(2));
    height = gelInfo.pocket.bounding_box(4);
    pocket_y = zeros(height+1,1);
    gelInfo.max_staple_migrate = 0 ;
    
    for i = 1:length(gelInfo.loop_indices)
    % starting with mono and then scaffold and ladder
    % (extrapolating staple migration in ladder and sacaffold we first need mono staple migration)
    
        indices = fliplr(gelInfo.loop_indices); 
        idx = indices(i);
        [spec,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        % Pocket fit
        pocket_y = pocket_y + spec.fullprofiles{pos}(pocket_position : pocket_position+height);
        
        % staples fits
        staple_fit = get_staple_fit(spec, pos, stapleHeight);
        gelInfo.species.(spec.type).staple.fits(pos, :) = staple_fit;
  
        % Ladder and Scaffold staple extrapolation from monomere staple migration
        try
            ends = gelInfo.species.mono.staple.migration_distance([1, end]);
            semi_ends = gelInfo.species.mono.staple.migration_distance([2, end-1]);
            mono_staple_migrate = semi_ends - ends;
        catch
        end
        
        if spec.type == "ladder"
            gelInfo.species.ladder.staple.migration_distance(pos) = ends(pos) - 2*mono_staple_migrate(pos);
            
        elseif spec.type == "scaffold"
            gelInfo.species.scaffold.staple.migration_distance(pos) = ends(pos) - mono_staple_migrate(pos);
            
        else
            gelInfo.species.mono.staple.migration_distance(pos) = staple_fit(2) - pocket_position;
        end
        
        
        % Species band fits
        [band_fit, fit_range] = get_band_fit(spec,  pos);
        gelInfo.species.(spec.type).fits(pos, :) = band_fit;
        gelInfo.species.(spec.type).fit_range(pos, :) = fit_range;
        gelInfo.species.(spec.type).band_positions(pos) = band_fit(2)';
        gelInfo.species.(spec.type).band_width(pos) = band_fit(3)';
        migration_distance = band_fit(2)' - pocket_position';
        
        % Migration distance
        gelInfo.species.(spec.type).migration_distance(pos) = migration_distance;
        
    end
    
    % Normalized Migration distance and band spread
    gelInfo = normalization(gelInfo);
  
    % Pocket fit
    pocket_x = double(pocket_position:pocket_position+height);
    pocket_fits =  coeffvalues(fit(pocket_x', pocket_y, 'gauss1'));
    gelInfo.pocket.fits = pocket_fits; 
    gelInfo.pocket.position = pocket_fits(2);
    gelInfo.pocket.width = pocket_fits(3);
    

    %% Display results and ask if ok
    close all
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot_band_fits(gelData, gelInfo)
    title(['Band positions with sigma ' num2str(gelInfo.sigma_integrate)])

    gelInfo.peaks_ok = strcmp(questdlg('Are the found peaks ok?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');            
    close all
        
end