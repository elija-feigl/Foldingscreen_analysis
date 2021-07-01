function gelInfo = integrate_species(gelInfo) 
% @ step2
% compute fractions of selected species per lane

    % integrate aggregates, monomers and smear
    n = length(gelInfo.loop_indices);

    
    mu = gelInfo.pocket.position;
    sig = gelInfo.pocket.width * gelInfo.sigma_integrate;

    for i = 1:n
        idx = gelInfo.loop_indices(i);
        [spec,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        if spec.type == "ladder" || spec.type == "scaffold"
            % integrate aggregates
            pocket_sums = 1.0;
            
            % integrate monomers
            monomer_sums = 0.0;
            
            % intergrate smear
            smear_sums = 0.0;
        else
            % integrate aggregates
            pocket_boundaries = [max(1,round(mu-sig)) round(mu+sig)];
            pocket_sums = sum(spec.fullprofiles{pos}(pocket_boundaries(1):pocket_boundaries(2) ));
            
            % integrate monomers
            mu_mono = spec.band_positions(pos);
            sig_mono = spec.band_width(pos) * gelInfo.sigma_integrate;
            monomer_boundaries = [max(1,round(mu_mono-sig_mono)) round(mu_mono+sig_mono)];
            monomer_sums = sum(spec.fullprofiles{pos}(monomer_boundaries(1):monomer_boundaries(2)));
            
            % intergrate smear
            smear_boundaries = [pocket_boundaries(2) monomer_boundaries(1)];
            smear_sums = sum(spec.fullprofiles{pos}(smear_boundaries(1):smear_boundaries(2)));
        end
            
        % Total bands
        total_band = smear_sums + monomer_sums + pocket_sums;
            
        % Fractions
        gelInfo.species.(spec.type).fraction_pocket(pos) = pocket_sums ./ total_band ;
        gelInfo.species.(spec.type).fraction_monomer(pos) = monomer_sums ./ total_band ;
        gelInfo.species.(spec.type).fraction_smear(pos) = smear_sums ./ total_band ;

    end
   
end

