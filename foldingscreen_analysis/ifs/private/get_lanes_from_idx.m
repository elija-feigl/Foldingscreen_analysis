    function [species,  pos] = get_lanes_from_idx(gelInfo, idx)
    % recieves the index from a linst of species indices (loop_indices) and returns the specie and position(in the lanes list) of
    % that index
    
        spec = gelInfo.species;
        if any(ismember(spec.ladder.indices, idx))
            species = spec.ladder;
            
        elseif any(ismember(spec.scaffold.indices, idx))
            species = spec.scaffold;
            
        else
            species = spec.mono;
        end
        
        pos = find(species.indices == idx); 
    end