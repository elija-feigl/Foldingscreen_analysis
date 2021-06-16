function [profileData, gelInfo] = select_species(profileData, gelData, gelInfo) %TODO: delete profileData
% @ step2
% select context specfic areas in gel and calculate bands with gauss fit

    mono_index = gelInfo.species.mono.indices;
    nLanes_mono =  length(mono_index);

    scaffold_index = gelInfo.species.scaffold.indices;
    nLanes_scaffold = length(scaffold_index);

    ladder_index = gelInfo.species.ladder.indices;
    nLanes_ladder = length(ladder_index);
    
    while ~gelInfo.peaks_ok
        
        lad_scaf_ok = false;
        mono_ok = false;
        staple_ok = false;
        pocket_ok = false;
        
        %% Pocket Selection
        plot_image_ui(gelData.images{1});
        hold on
        
        while ~pocket_ok
            % select area for pockets
            title('Select pockets area')
            rect_pocket = drawrectangle('Label','Pocket','Color',[1 0 0]);
            gelInfo.pocket.bounding_box = int32(rect_pocket.Position);
            pocket_ok = strcmp(questdlg('Are you ok with the selected Pocket?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');
            if ~pocket_ok
                delete(rect_pocket)
            end
        end
        
        %% Ladder and Scaffold selection
        while ~lad_scaf_ok
            lad_scaf_indices = sort([ladder_index  scaffold_index]);
            ladder_counter = 0;
            scaffold_counter = 0;
            
            for i=1:length(lad_scaf_indices)
                if ismember(lad_scaf_indices(i), ladder_index)
                    ladder_counter = ladder_counter + 1; % indentation on the gel image
                    title(sprintf('Select Ladder %.f', ladder_counter))
                    pos_ladder = drawpoint('Label',sprintf('L %.f', ladder_counter),'Color','y');
                    %NOTE: we only take the y component and discard the x component for now
                    gelInfo.species.ladder.positions(ladder_counter) = int32(pos_ladder.Position(2))';
                    
                    
                elseif ismember(lad_scaf_indices(i), scaffold_index)
                    scaffold_counter = scaffold_counter + 1; % indentation on the gel image
                    title(sprintf('Select Scaffold %.f', scaffold_counter))
                    pos_scaf = drawpoint('Label',sprintf('S %.f', scaffold_counter),'Color',[0 1 0]);
                    %NOTE: we only take the y component and discard the x component for now
                    gelInfo.species.scaffold.positions(scaffold_counter) = int32(pos_scaf.Position(2))';
                    

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
            
            pos_mono = line_mono.Position(:, 2);
            if length(pos_mono) == 2
                pos_mono = linspace(pos_mono(1), pos_mono(2), nLanes_mono);
            else
                pos_mono = int32(pos_mono);
            end
            %NOTE: we only take the y component and discard the x component for now
            gelInfo.species.mono.positions = pos_mono(1:nLanes_mono);

        end
        
        %% Staple Selection
        
        while ~staple_ok
            title('Select staple line (double click to place last)')
            line_staple = drawpolyline('Label','Staple','Color',[0.5 0.5 0.5]);

            staple_ok = strcmp(questdlg('Are you ok with the selected points?','Peaks found?' ,'No','Yes', 'Yes'),'Yes');
            if ~staple_ok
                delete(line_staple)
                continue
            end

            pos_staple = line_staple.Position;
            
            %NOTE: we only take the y component and discard the x component for now
            if length(pos_staple) == 2
                selectedStaple_y = linspace(pos_staple(1,2), pos_staple(2,2), nLanes_mono);
                gelInfo.species.mono.staple.positions = selectedStaple_y';

            elseif length(pos_staple) == 1
                selectedStaple_y = linspace(pos_staple(1,2), pos_staple(1,2), nLanes_mono);
                gelInfo.species.mono.staple.positions = selectedStaple_y';

            else
                pos_staple = int32(pos_staple);
                gelInfo.species.mono.staple.positions = pos_staple(1:nLanes_mono, :);
            end
            gelInfo.species.ladder.staple.positions = zeros(nLanes_ladder);
            gelInfo.species.scaffold.staple.positions = zeros(nLanes_scaffold);

        end
        
        %TODO: -low- add listener to update positions
        title('double click pocket-area to finish)')
        wait(rect_pocket);
        close all
        
        [profileData, gelInfo] = species_fits(profileData, gelInfo, gelData);
    end
end  