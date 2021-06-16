        
%% Calculate fits
function [profileData, gelInfo] = species_fits(profileData, gelInfo, gelData) %TODO: delete profileData

    %% Get band fits
    function fits = get_band_fit(species, pos)
        xpos = species.positions(pos) - fix(species.height/2);
        fits = zeros(1, 3);
        
        upper = xpos + species.height;
        val = species.fullprofiles{pos};
        
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
    
    %TODO: delete after adapting part 3
    mono_index = gelInfo.species.mono.indices;
    nLanes_mono =  length(mono_index); 
    
    %TODO: delete after adapting part 3
    scaffold_index = gelInfo.species.scaffold.indices;
    nLanes_scaffold = length(scaffold_index); 
    
    %TODO: delete after adapting part 3
    ladder_index = gelInfo.species.ladder.indices;
    nLanes_ladder = length(ladder_index);  
    
    %TODO: delete fits, after adapting part 3
    gelInfo.species.ladder.fits = zeros(nLanes_ladder, 3);
    gelInfo.species.scaffold.fits = zeros(nLanes_scaffold, 3);
    gelInfo.species.mono.fits = zeros(nLanes_mono, 3);

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
    
    for i = 1:length(gelInfo.loop_indices)
    % starting with mono and then scaffold and ladder
    % (extrapolating staple migration in ladder and sacaffold we first need mono staple migration)
    
        indices = fliplr(gelInfo.loop_indices); 
        idx = indices(i);
        [species,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        % Pocket fit
        pocket_y = pocket_y + species.fullprofiles{pos}(pocket_position : pocket_position+height);
        
        % Species fit
        band_fit = get_band_fit(species,  pos);
        gelInfo.species.(species.type).fits(pos, :) = band_fit; %TODO: delete after adapting part 3
        gelInfo.species.(species.type).positions(pos) = band_fit(2)';
        gelInfo.species.(species.type).band_width(pos) = band_fit(3)';
        gelInfo.species.(species.type).migration_distance(pos) = band_fit(2)' - pocket_position';
        
        % For staples
        staple_fit = get_staple_fit(species, pos, stapleHeight);
        gelInfo.species.(species.type).staple.fits(pos, :) = staple_fit; %TODO: delete after adapting part 3
  
        
        % Ladder and Scaffold staple extrapolation from monomere staple migration
        try
            ends = gelInfo.species.mono.staple.migration_distance([1, end]);
            semi_ends = gelInfo.species.mono.staple.migration_distance([2, end-1]);
            mono_staple_migrate = semi_ends - ends;
        catch
        end
        
        if species.type == "ladder"
            gelInfo.species.ladder.staple.migration_distance(pos) = ends(pos) - 2*mono_staple_migrate(pos);
            
        elseif species.type == "scaffold"
            gelInfo.species.scaffold.staple.migration_distance(pos) = ends(pos) - mono_staple_migrate(pos);
            
        else
            gelInfo.species.mono.staple.migration_distance(pos) = staple_fit(2) - pocket_position;
        end
        
    end
  
    % Pocket fit
    pocket_x = double(pocket_position:pocket_position+height);
    pocket_fits =  coeffvalues(fit(pocket_x', pocket_y, 'gauss1'));
    gelInfo.pocket.fits = pocket_fits; 
    gelInfo.pocket.position = pocket_fits(2);
    gelInfo.pocket.width = pocket_fits(3);
    
    %% get band data

    % TODO: delete after adapting part 3
    mono_fits = zeros(nLanes_ladder+nLanes_mono+nLanes_scaffold, 3);
    
    if ~isempty(gelInfo.species.ladder.indices)
        mono_fits(gelInfo.species.ladder.indices, :) = gelInfo.species.ladder.fits;
    end

    if ~isempty(gelInfo.species.scaffold.indices)   
        mono_fits(gelInfo.species.scaffold.indices, :) = gelInfo.species.scaffold.fits;
    end

    mono_fits(gelInfo.species.mono.indices, :) = gelInfo.species.mono.fits;


    staples_fits = zeros(nLanes_ladder + nLanes_mono+ nLanes_scaffold, 3);
    staples_fits(gelInfo.species.mono.indices, :) = gelInfo.species.mono.staple.fits;

    %% add to profiles structure
    % NOTE: output necessary for legacy python code(ifs database code) 
    % TODO: delete after adapting part 3
    profileData.monomerFits = mono_fits; % ladder, scaffold and monomers fits
    profileData.stapleFits = staples_fits;  

    %% Display results and ask if ok
    close all
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot_band_fits(gelData, gelInfo)
    title(['Band positions with sigma ' num2str(gelInfo.sigma_integrate)])

    gelInfo.peaks_ok = strcmp(questdlg('Are the found peaks ok?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');            
    close all
        
end