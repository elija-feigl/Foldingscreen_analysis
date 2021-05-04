function [profileData, gelInfo] = select_species(profileData, gelData, gelInfo)
% @ step2
% select context specfic areas in gel and calculate bands with gauss fit
% TODO optimize selection order for fast & convenient workflow

    function fits = calculate_fit(profileData, species, height, xpos_pocket)
        if species.type ~= "none"
            index = species.index;
            xpos = species.positions(: ,2) - fix(height/2);
            fits = zeros(length(index), 3);
            
            for j=1:length(xpos)
                upper = xpos(j)+height;
                if species.type == "staple"
                    if upper > length(profileData.fullProfiles{index(j)})
                        upper = length(profileData.fullProfiles{index(j)});
                    end
                end
                val = profileData.fullProfiles{index(j)};
                if xpos(j) < xpos_pocket
                    fits(j,1) = 0.0;
                    fits(j,2) = xpos_pocket;
                    fits(j,3) = 0.0;
                else
                    y = val(xpos(j):upper);
                    x = double(xpos(j):upper);
                    options = fitoptions('gauss1', 'Lower', [0 xpos(j) 0], 'Upper', [Inf upper Inf]);

                    fits(j,:) = coeffvalues(fit(x', y, 'gauss1', options));
                end

            
            end
        else
            fits = zeros(2,3);
            fits(:,2) = xpos_pocket;
        end
    end
    
    n = length(profileData.fullProfiles);    
    peaks_ok = false;

    % TODO hardcoding these might not be ideal, as pixelsize can depend
    %       on the scanner used. find better solution!
    ladderHeight = 50;
    scaffoldHeight = 50;
    bandHeight = 50;
    stapleHeight = 200; 
    
    mono_index = gelInfo.species.mono.index;
    nLanes_mono =  length(mono_index);
    
    scaffold_index = gelInfo.species.scaffold.index;
    nLanes_scaffold = length(scaffold_index);
    
    ladder_index = gelInfo.species.ladder.index;
    nLanes_ladder = length(ladder_index);
    
    if length(ladder_index) > 1
        has_ladder = true;
    else
        has_ladder = false;
    end
    
    while ~peaks_ok
        
        lad_scaf_ok = false;
        mono_ok = false;
        staple_ok = false;
        pocket_ok = false;
        
        %% Pocket Selection
        plot_image_ui(gelData.images{1});
        hold on
        
        % select area for pockets
        while ~pocket_ok
            title('Select pockets area')
            rect_pocket = drawrectangle('Label','Pocket','Color',[1 0 0]);
            gelInfo.species.pocket.positions = int32(rect_pocket.Position);
            pocket_ok = strcmp(questdlg('Are you ok with the selected Pocket?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');
            if ~pocket_ok
                delete(rect_pocket)
            end
        end
        
        %% Ladder and Scaffold selection
        while ~lad_scaf_ok
            lad_scaf_indeces = sort([ladder_index  scaffold_index]);
            ladder_counter = 0;
            scaffold_counter = 0;
            for i=1:length(lad_scaf_indeces)
                if ismember(lad_scaf_indeces(i), ladder_index)
                    ladder_counter = ladder_counter + 1; % indentation on the gel image
                    title(sprintf('Select Ladder %.f', ladder_counter))
                    pos_ladder = drawpoint('Label',sprintf('L %.f', ladder_counter),'Color','y');
                    gelInfo.species.ladder.positions(ladder_counter, :) = int32(pos_ladder.Position);
                    
                elseif ismember(lad_scaf_indeces(i), scaffold_index)
                    scaffold_counter = scaffold_counter + 1; % indentation on the gel image
                    title(sprintf('Select Scaffold %.f', scaffold_counter))
                    pos_scaf = drawpoint('Label',sprintf('S %.f', scaffold_counter),'Color',[0 1 0]);
                    gelInfo.species.scaffold.positions(scaffold_counter, :) = int32(pos_scaf.Position);
                else
                    disp("Design has no Ladder nor Scaffold")
                end
            end
            lad_scaf_ok = strcmp(questdlg('Are you ok with the selected points?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');
            if ~lad_scaf_ok
                delete(pos_scaf)
                delete(pos_ladder)
                continue
            end
        end

        %% Monomer Selection
        while ~mono_ok
            title('Select monomer Bands. Double click to place last and pick faulty lanes above pocket to ignore')
            line_mono = drawpolyline('Label','Monomer','Color',[0 0 1]);
            mono_ok = strcmp(questdlg('Are you ok with the selected points?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');
            if ~mono_ok
                delete(line_mono)
                continue
            end
            
            pos_mono = line_mono.Position;
            if length(pos_mono) == 2
                selectedMono_y = linspace(pos_mono(1,2), pos_mono(2,2), nLanes_mono);
                pos_mono = [ones(1,nLanes_mono,'uint32')' selectedMono_y'];
            else
                pos_mono = int32(pos_mono);
            end

            gelInfo.species.mono.positions = pos_mono(1:nLanes_mono, :);
        end
        
        %% Staple Selection
        
        gelInfo.species.staples.index = gelInfo.species.mono.index;
        gelInfo.species.staples.type = "staple";
        
        while ~staple_ok
            title('Select staple line (double click to place last)')
            line_staple = drawpolyline('Label','Staple','Color',[0.5 0.5 0.5]);

            staple_ok = strcmp(questdlg('Are you ok with the selected points?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');
            if ~staple_ok
                delete(line_staple)
                continue
            end

            pos_staple = line_staple.Position;
            
            if length(pos_staple) == 2
                selectedStaple_y = linspace(pos_staple(1,2), pos_staple(2,2), nLanes_mono);
                gelInfo.species.staples.positions = [ones(1,nLanes_mono,'uint32')' selectedStaple_y'];

            elseif length(pos_staple) == 1
                selectedStaple_y = linspace(pos_staple(1,2), pos_staple(1,2), nLanes_mono);
                gelInfo.species.staples.positions = [ones(1,nLanes_mono,'uint32')' selectedStaple_y'];

            else
                pos_staple = int32(pos_staple);
                gelInfo.species.staples.positions = pos_staple(1:nLanes_mono, :);
            end
        end
        
        %TODO: -low- add listener to update positions
        title('double click pocket-area to finish)')
        wait(rect_pocket);
        close all
        %% Calculate fits
        % compute pocket sum profiles and fit it with gaussian.
        % position and width of pocket is always the same -> sum
        
        xpos_pocket = gelInfo.species.pocket.positions(2);
        height = gelInfo.species.pocket.positions(4);
        y = zeros(height+1,1);
        for i=1:n
            y = y + profileData.fullProfiles{i}(xpos_pocket:xpos_pocket+height);
        end
        x = double(xpos_pocket:xpos_pocket+height);
        gelInfo.species.pocket.fits =  coeffvalues(fit(x', y, 'gauss1'));
        
        gelInfo.species.ladder.fits = calculate_fit(profileData, gelInfo.species.ladder, ladderHeight, xpos_pocket); 
        gelInfo.species.scaffold.fits = calculate_fit(profileData, gelInfo.species.scaffold, scaffoldHeight, xpos_pocket);
        gelInfo.species.mono.fits = calculate_fit(profileData, gelInfo.species.mono, bandHeight, xpos_pocket); 
        gelInfo.species.staples.fits = calculate_fit(profileData, gelInfo.species.staples, stapleHeight, xpos_pocket);
        
        %% get band data
        
        % TODO: change this after adapting the file: plot_band_fits to the
        % new changes
        mono_fits = zeros(nLanes_ladder+nLanes_mono+nLanes_scaffold, 3);
        if gelInfo.species.ladder.type ~= "none"
            mono_fits(gelInfo.species.ladder.index, :) = gelInfo.species.ladder.fits;
        end
        if gelInfo.species.scaffold.type ~= "none"   
            mono_fits(gelInfo.species.scaffold.index, :) = gelInfo.species.scaffold.fits;
        end
            
        mono_fits(gelInfo.species.mono.index, :) = gelInfo.species.mono.fits;
        
        
        staples_fits = zeros(nLanes_ladder + nLanes_mono+ nLanes_scaffold, 3);
        staples_fits(gelInfo.species.mono.index, :) = gelInfo.species.staples.fits;
       
        %% add to profiles structure
        % NOTE: output necessary for legacy python code(ifs database code) 
        
        profileData.aggregateFit = gelInfo.species.pocket.fits;
        profileData.aggregateSelectedArea = gelInfo.species.pocket.positions;
        profileData.monomerFits = mono_fits;
        profileData.monomerSelectedArea = gelInfo.species.mono.positions;
        profileData.stapleLine = gelInfo.species.staples.positions;
        profileData.stapleFits = staples_fits;
        profileData.has_ladder = has_ladder;
    
        %% Display results and ask if ok
        close all
        figure('units','normalized','outerposition',[0 0 1 1]);
        plot_band_fits(gelData, profileData, gelInfo)
        title(['Band positions with sigma *' num2str(profileData.sigma_integrate)])
        
        peaks_ok = strcmp(questdlg('Are the found peaks ok?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');            
        close all
    end
end  