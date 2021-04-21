function profileData = select_species(profileData,gelData, gelInfo)
% @ step2
% select context specfic areas in gel and calculate bands with gauss fit
% TODO optimize selection order for fast & convenient workflow
    n = length(profileData.fullProfiles);    
    peaks_ok= false;
    while ~peaks_ok
        %% select species
        plot_image_ui(gelData.images{1});
        % TODO use gelInfo to figure out how many folding lanes to select
        % TODO use gelInfo to figure out how many scaffolds to select
        % TODO use gelInfo to figure out how many ladders to select
        nLanes_mono = 16;
        nLanes_scaffold = 2;
        nLanes_ladder = 2;
        % select area for pockets
        title('Select pockets area')
        rect_pocket = drawrectangle('Label','Pocket','Color',[1 0 0]);
        selectedPocketArea = int32(rect_pocket.Position);
        % select line for monomer band
        title('Select monomer line (double click to place last). pick faulty lanes above pocket to ignore')
        line_mono = drawpolyline('Label','Monomer','Color',[0 0 1]);
        pos_mono = line_mono.Position;
        if length(pos_mono) == 2
            selectedMono_y = linspace(pos_mono(1,2), pos_mono(2,2), nLanes_mono);
            selectedMono = [ones(1,nLanes_mono,'uint32')' selectedMono_y'];
        else
            selectedMono = int32(pos_mono);
        end   
        % select scaffold bands
        % TODO use nLanes_scaffold to figure out how many scaff to select
        title('Select Scaffold L')
        point_scaffL = drawpoint('Label','sL','Color',[0 1 0]);
        selectedScaffold_L = int32(point_scaffL.Position);
        title('Select Scafffold R')
        point_scaffR = drawpoint('Label','sR','Color',[0 1 0]);
        selectedScaffold_R = int32(point_scaffR.Position);
        % select ladder bands
        % TODO instead of the prompt use nLanes_ladder
        has_ladder = strcmp(questdlg('Ladder?','Does it include a ladder?' ,'No','Yes', 'Yes'),'Yes');
        if has_ladder    
            title('Select ladder L')
            point_ladderL = drawpoint('Label','lL','Color',[0 1 0]);
            selectedLadder_L = int32(point_ladderL.Position);
            title('Select ladder R')
            point_ladderR = drawpoint('Label','lR','Color',[0 1 0]);
            selectedLadder_R = int32(point_ladderR.Position);
        end
        % select line for staple band
        title('Select staple line (double click to place last)')
        line_staple = drawpolyline('Label','Staple','Color',[0.5 0.5 0.5]);
        pos_staple = line_staple.Position;
        if length(pos_staple) == 2
            selectedStaple_y = linspace(pos_staple(1,2), pos_staple(2,2), nLanes_mono);
            selectedStaple = [ones(1,nLanes_mono,'uint32')' selectedStaple_y'];
            
        elseif length(pos_staple) == 1
            selectedStaple_y = linspace(pos_staple(1,2), pos_staple(1,2), nLanes_mono);
            selectedStaple = [ones(1,nLanes_mono,'uint32')' selectedStaple_y'];
            
        else
            selectedStaple = int32(pos_staple);
        end
        
        %NOTE: positions that have been moved around do not update
        %TODO: -low- add listener to update positions
        wait(rect_pocket);
        title('double click pocket-area to finish)')
        close all
        
        %% get band data
        % compute pocket sum profiles and fit it with gaussian.
        %       position and width of pocket is always the same -> sum
        xpos_pocket = selectedPocketArea(2);
        height = selectedPocketArea(4);
        y = zeros(height+1,1);
        for i=1:n
            y = y + profileData.fullProfiles{i}(xpos_pocket:xpos_pocket+height);
        end
        x = double(xpos_pocket:xpos_pocket+height);
        pocket_fit =  coeffvalues(fit(x', y, 'gauss1'));

        % compute monomer profiles and fit it with gaussian
        % TODO hardcoding these might not be ideal, as pixelsize can depend
        %       on the scanner used. find better solution!
        ladderHeight = 50;
        scaffoldHeight = 50;
        bandHeight = 50;
        stapleHeight = 200; 
        mono_fits = zeros(length(profileData.fullProfiles),3);
        staple_fits = zeros(length(profileData.fullProfiles),3);
        
        for i=1:n     
            if has_ladder
                if i == 1
                    height = ladderHeight;
                    xpos = selectedLadder_L(2) - fix(height/2);
                elseif i == 2
                    height = scaffoldHeight;
                    xpos = selectedScaffold_L(2) - fix(height/2);
                elseif i == (n-1)
                    height = scaffoldHeight;
                    xpos = selectedScaffold_R(2) - fix(height/2);
                elseif i == n
                    height = ladderHeight;
                    xpos = selectedLadder_R(2) - fix(height/2);
                else
                    height = bandHeight;
                    xpos =  selectedMono(i-2,2) - fix(height/2);
                    xpos_staple =  selectedStaple(i-2,2) - fix(height/2);
                end
            else
                if i == 1
                    height = scaffoldHeight;
                    xpos = selectedScaffold_L(2) - fix(height/2);
                elseif i == n
                    height = scaffoldHeight;
                    xpos = selectedScaffold_R(2) - fix(height/2);
                else
                    height = bandHeight;
                    xpos =  selectedMono(i-1,2) - fix(height/2);
                    xpos_staple =  selectedStaple(i-1,2) - fix(height/2);
                end
            end
            y = profileData.fullProfiles{i}(xpos:xpos+height);
            x = double(xpos:xpos+height);
            options = fitoptions('gauss1', 'Lower', [0 xpos 0], 'Upper', [Inf (xpos+height) Inf]);
            mono_fits(i,:) = coeffvalues(fit(x', y, 'gauss1', options)); 
            
            if i > has_ladder+1 && i < n-has_ladder
                upper = xpos_staple+stapleHeight;
                if upper > length(profileData.fullProfiles{i})
                    upper = length(profileData.fullProfiles{i});
                end

                y = profileData.fullProfiles{i}(xpos_staple:upper);
                x = double(xpos_staple:upper);
                options = fitoptions('gauss1', 'Lower', [0 xpos_staple 0], 'Upper', [Inf upper Inf]);
                staple_fits(i,:) = coeffvalues(fit(x', y, 'gauss1', options)); 
            
            end
            % ignore lanes picked above pocket
            if xpos < xpos_pocket
                mono_fits(i,1) = 0.0;
                mono_fits(i,2) = xpos_pocket;
                mono_fits(i,3) = 0.0;
            end
        end    

        % add to profiles structure
        profileData.aggregateFit = pocket_fit;
        profileData.aggregateSelectedArea = selectedPocketArea;
        profileData.monomerFits = mono_fits;
        profileData.monomerSelectedArea = selectedMono;
        profileData.stapleLine = selectedStaple;
        profileData.stapleFits = staple_fits;
        profileData.has_ladder = has_ladder;
    
        %% Display results and ask if ok
        close all
        figure('units','normalized','outerposition',[0 0 1 1]);
        plot_band_fits(gelData, profileData)
        title(['Band positions with sigma *' num2str(profileData.sigma_integrate)])
        
        peaks_ok = strcmp(questdlg('Are the found peaks ok?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');
        close all
    end
end