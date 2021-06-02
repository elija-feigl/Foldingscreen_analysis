function [profileData] = integrate_species(profileData, gelInfo)
% @ step2
% compute fractions of selected species per lane

    % integrate aggregates, monomers and smear
    n = length(gelInfo.loop_indices);
    
    mu = gelInfo.pocket.fits(2);
    sig = gelInfo.pocket.fits(3) * gelInfo.sigma_integrate;
    
    pocket_sums = zeros(n, 1);
    pocket_boundaries = zeros(n,2);
    
    monomer_sums = zeros(n,1);
    monomer_boundaries = zeros(n,2);
    
    smear_sums = zeros(n,1);
    smear_boundaries = zeros(n,2);

    for i = 1:n
        idx = gelInfo.loop_indices(i);
        [species,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        if species.type == "ladder" || species.type == "scaffold"
            % integrate aggregates
            pocket_sums(idx) = 1.0;
            
            % integrate monomers
            monomer_sums(idx) = 0.0;
            
            % intergrate smear
            smear_sums(i) = 0.0;
        else
            % integrate aggregates
            pocket_boundaries(idx,:) = [max(1,round(mu-sig)) round(mu+sig)];
            pocket_sums(idx) = sum(species.fullprofiles{pos}(pocket_boundaries(idx,1):pocket_boundaries(idx,2) ));
            
            % integrate monomers
            mu_mono = species.fits(pos,2);
            sig_mono = species.fits(pos,3) * gelInfo.sigma_integrate;
            monomer_boundaries(idx,:) = [max(1,round(mu_mono-sig_mono)) round(mu_mono+sig_mono)];
            monomer_sums(idx) = sum(species.fullprofiles{pos}(monomer_boundaries(idx,1):monomer_boundaries(idx,2)));
            
            % intergrate smear
            smear_boundaries(idx,:) = [pocket_boundaries(idx,2) monomer_boundaries(idx,1)];
            smear_sums(idx) = sum(species.fullprofiles{pos}(smear_boundaries(idx,1):smear_boundaries(idx,2)));
        end
    end
    
    profileData.monomerBoundaries = monomer_boundaries; %TODO: not used
    profileData.pocketBoundaries = pocket_boundaries; %TODO: not used
    profileData.monomerTotal = monomer_sums;
    profileData.pocketTotal = pocket_sums;
    profileData.smearTotal = smear_sums;
    profileData.smearBoundaries = smear_boundaries; %TODO: not used
    
end

