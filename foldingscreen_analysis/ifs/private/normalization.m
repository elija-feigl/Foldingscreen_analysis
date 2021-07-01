function gelInfo = normalization(gelInfo)
% Compute normalized band_width, mono_spread, migrations
    
    keySet = {'1317' , '2873' , '4536' , '7249' , '7560' , '7704' , '8064' , '9072' , 'CS11' , 'CS12' , 'CS13' , 'CS15' , 'CS16' , 'CS17'};
    valueSet = {0.949670, 0.757241, 0.617800, 0.521126, 0.524794, 0.573690, 0.504229, 0.412140, 0.473944, 0.435950, 0.537250, 0.518357, 0.490852, 0.414030};
    scaffold_refactor = containers.Map (keySet, valueSet);


    n = length(gelInfo.loop_indices);
   
    %% calculate overall nomalization factor
    Norm_factor = 0;
    max_staple_migrate = max(gelInfo.species.mono.staple.migration_distance);
    
    for i=1:n
        idx = gelInfo.loop_indices(i);
        [spec,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        %NOTE: sqrt of stapleNorm to reduce its effect on the normalisation.Full effect seems to overcorrect migration
        %distances (theoretical explanation missing!)
        gelNorm = sqrt(gelInfo.species.(spec.type).staple.migration_distance(pos) ./ max_staple_migrate);
        
        if spec.type == "ladder" && Norm_factor == 0
            Norm_factor = gelInfo.species.ladder.migration_distance(pos) * gelNorm;
             
        elseif spec.type == "scaffold" && Norm_factor == 0
            Norm_factor = gelInfo.species.scaffold.migration_distance(pos) ./ scaffold_refactor(gelInfo.scaffold_type) * gelNorm;
            
        elseif spec.type == "mono" && Norm_factor == 0
            disp('No Ladder nor Scaffold for Normalization') %TODO: exit code
        end
        
        gelInfo.species.(spec.type).band_spread(pos) = gelInfo.species.(spec.type).band_width(pos) ./ Norm_factor ./ gelNorm;
        gelInfo.species.(spec.type).migration_distance(pos) = gelInfo.species.(spec.type).migration_distance(pos) ./ Norm_factor ./ gelNorm ;
  
    end
  
end
