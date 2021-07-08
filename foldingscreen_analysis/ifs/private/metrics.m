function gelInfo = metrics(gelInfo)    
    %%New Metric method 
   spreadNormfactor = 15.0;
   mono_migrate_best = max(gelInfo.species.mono.migration_distance);
   mono_spread_best = max(gelInfo.species.mono.band_spread);
   
   n = length(gelInfo.loop_indices);
   
   for i = 1:n
        idx = gelInfo.loop_indices(i);
        [spec,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        if spec.type == "ladder" || spec.type == "scaffold"
            gelInfo.species.(spec.type).folding_quality_metric(pos) = 0.0;
        else
            %NOTE: assuming band spread bigger than one means double peak or a bad band
            if spec.band_spread(pos) > 1.0 / spreadNormfactor
                folding_quality_metric = 0;
            else
                folding_quality_metric = spec.fraction_monomer(pos) .* (1.0 - spec.band_spread(pos) .* spreadNormfactor);
            end
            gelInfo.species.(spec.type).folding_quality_metric(pos) = folding_quality_metric;            
        end
        
        % compute relative migration distance
        gelInfo.species.(spec.type).rel_band_migrate(pos) = gelInfo.species.(spec.type).migration_distance(pos)./mono_migrate_best;
        gelInfo.species.(spec.type).rel_band_spread(pos) = gelInfo.species.(spec.type).band_spread(pos)./mono_spread_best;
        
   end
   
   %TODO: check the exact conditions to be correct
   % compute relative migration distance
   if length(gelInfo.species.ladder.indices) >= 2
        ladder  = gelInfo.species.ladder;
        gelInfo.species.ladder.migrate_error = abs(ladder.rel_band_migrate(1) - ladder.rel_band_migrate(end));
        gelInfo.species.ladder.spread_error = abs(ladder.rel_band_spread(1) - ladder.rel_band_spread(end));
   else
       gelInfo.species.ladder.migrate_error = 0.1 ; 
       gelInfo.species.ladder.spread_error = 0.1 ; 
   end
   
   %NOTE: scaffold error less accurate - > factor 2
   if length(gelInfo.species.ladder.indices) < 2 && length(gelInfo.species.scaffold.indices) >= 2
        scaffold  = gelInfo.species.scaffold;
        gelInfo.species.scaffold.migrate_error = 2.0 * abs(scaffold.rel_band_migrate(1) - scaffold.rel_band_migrate(end));
        gelInfo.species.scaffold.spread_error = 2.0 * abs(scaffold.rel_band_spread(1) - scaffold.rel_band_spread(end));
   else
       gelInfo.species.scaffold.migrate_error = 0.1;
       gelInfo.species.scaffold.spread_error = 0.1 ; 
   end
   

    cutoff = 0.75; 
    if ~isempty(gelInfo.species.ladder.indices)
        tolerance = 1.0 * gelInfo.species.ladder.migrate_error;
    else
        tolerance = 1.0 * gelInfo.species.scaffold.migrate_error;
    end
    
    %NOTE: 0.03 shown to be reasonable lower limit
    if tolerance < 0.03
        tolerance = 0.03;
    end
    
    %% Folding quality metric over migration
    for i = 1:n
        idx = gelInfo.loop_indices(i);
        [spec,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        rel_band_migrate = gelInfo.species.(spec.type).rel_band_migrate(pos);
        folding_quality_metric = gelInfo.species.(spec.type).folding_quality_metric(pos);
        
        if rel_band_migrate < cutoff
            gelInfo.species.(spec.type).folding_quality_metric_migrate(pos) = 0.0;
        elseif abs(rel_band_migrate - 1.0) < tolerance
            gelInfo.species.(spec.type).folding_quality_metric_migrate(pos) = folding_quality_metric;
        else
            %NOTE: ^2 increasses the harsheness of the migration penalty to balance the migration quality towards quality 
            %(no explanation theoretical explanation given)
            migrate_penalty = (rel_band_migrate - cutoff) / (1.0 - cutoff);
            gelInfo.species.(spec.type).folding_quality_metric_migrate(pos) = migrate_penalty^2 .* folding_quality_metric; 
        end
    end
end