        
%% Calculate fits
function [profileData, gelInfo, peaks_ok] = species_fits(profileData, gelInfo, gelData)

    %% Get band fits
    function fits = get_band_fit(species, pos)
        xpos = species.positions(pos) - fix(species.height/2);
        fits = zeros(1, 3);
        
        upper = xpos + species.height;
        val = species.fullprofiles{pos};
        
        % pocket_position = gelInfo.pocket.positions(2);
        
        if xpos < pocket_position
            fits(1) = 0.0;
            fits(2) = pocket_position;
            fits(3) = 0.0;
        else
            y = val(xpos:upper);
            x = double(xpos:upper);
            options = fitoptions('gauss1', 'Lower', [0 xpos 0], 'Upper', [Inf upper Inf]);

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
            xpos = species.positions(pos) - fix(height/2);
            upper = xpos + height;
            val = species.fullprofiles{pos};   

            if upper > length(species.fullprofiles{pos})
                upper = length(species.fullprofiles{pos});
            end

            y = val(xpos:upper);
            x = double(xpos:upper);
            options = fitoptions('gauss1', 'Lower', [0 xpos 0], 'Upper', [Inf upper Inf]);

            fits = coeffvalues(fit(x', y, 'gauss1', options));
        end
    end

    mono_index = gelInfo.species.mono.indices;
    nLanes_mono =  length(mono_index);

    scaffold_index = gelInfo.species.scaffold.indices;
    nLanes_scaffold = length(scaffold_index);

    ladder_index = gelInfo.species.ladder.indices;
    nLanes_ladder = length(ladder_index);

    gelInfo.species.ladder.fits = zeros(nLanes_ladder, 3);
    gelInfo.species.scaffold.fits = zeros(nLanes_scaffold, 3);
    gelInfo.species.mono.fits = zeros(nLanes_mono, 3);

    % TODO hardcoding these might not be ideal, as pixelsize can depend
    %       on the scanner used. find better solution!
    gelInfo.species.ladder.height = 50;
    gelInfo.species.scaffold.height = 50;
    gelInfo.species.mono.height = 50;
    stapleHeight = 200; 

    %% Pocket, Ladder, Scaffold, Monomers and staple fits
    
    % compute pocket sum profiles and fit it with gaussian.
    % position and width of pocket is always the same -> sum
    
    pocket_position = gelInfo.pocket.positions(2);
    height = gelInfo.pocket.positions(4);
    pocket_y = zeros(height+1,1);
    
    for i = 1:length(gelInfo.loop_indices)
        idx = gelInfo.loop_indices(i);
        [species,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        % Pocket fit
        pocket_y = pocket_y + species.fullprofiles{pos}(pocket_position:pocket_position+height);
        
        % Species fit
        gelInfo.species.(species.type).fits(pos, :) = get_band_fit(species,  pos);
        gelInfo.species.(species.type).staple.fits(pos, :) = get_staple_fit(species, pos, stapleHeight);
    end
    
    % Pocket fit
    pocket_x = double(pocket_position:pocket_position+height);
    gelInfo.pocket.fits =  coeffvalues(fit(pocket_x', pocket_y, 'gauss1'));
    
    %% get band data
    
    % TODO: change this after adapting the file: plot_band_fits to the
    % new changes
    mono_fits = zeros(nLanes_ladder+nLanes_mono+nLanes_scaffold, 3);
    
    if isempty(gelInfo.species.ladder.indices)
        mono_fits(gelInfo.species.ladder.indices, :) = gelInfo.species.ladder.fits;
    end

    if isempty(gelInfo.species.scaffold.indices)   
        mono_fits(gelInfo.species.scaffold.indices, :) = gelInfo.species.scaffold.fits;
    end

    mono_fits(gelInfo.species.mono.indices, :) = gelInfo.species.mono.fits;


    staples_fits = zeros(nLanes_ladder + nLanes_mono+ nLanes_scaffold, 3);
    staples_fits(gelInfo.species.mono.indices, :) = gelInfo.species.mono.staple.fits;

    %% add to profiles structure
    % NOTE: output necessary for legacy python code(ifs database code) 


    profileData.aggregateFit = gelInfo.pocket.fits;
    profileData.aggregateSelectedArea = gelInfo.pocket.positions;
    profileData.monomerFits = mono_fits; % ladder, scaffold and monomers fits
    profileData.monomerSelectedArea = gelInfo.species.mono.positions;
    profileData.stapleLine = gelInfo.species.mono.staple.positions;
    profileData.stapleFits = staples_fits; 
    

    %% Display results and ask if ok
    close all
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot_band_fits(gelData, gelInfo)
    title(['Band positions with sigma ' num2str(gelInfo.sigma_integrate)])

    peaks_ok = strcmp(questdlg('Are the found peaks ok?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');            
    close all
        
end