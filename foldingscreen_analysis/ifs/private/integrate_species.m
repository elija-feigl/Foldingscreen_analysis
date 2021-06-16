function [profileData, gelInfo] = integrate_species(profileData, gelInfo) %TODO: delete profileData
% @ step2
% compute fractions of selected species per lane

    % integrate aggregates, monomers and smear
    n = length(gelInfo.loop_indices);
    
    %TODO: delete after adapting part 3
    profileData.monomerTotal = zeros(n, 1);
    profileData.pocketTotal = zeros(n, 1);
    profileData.smearTotal = zeros(n, 1);
    
    
    mu = gelInfo.pocket.position;
    sig = gelInfo.pocket.width * gelInfo.sigma_integrate;

    for i = 1:n
        idx = gelInfo.loop_indices(i);
        [species,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        if species.type == "ladder" || species.type == "scaffold"
            % integrate aggregates
            gelInfo.species.(species.type).pocket_sums(pos) = 1.0;
            profileData.pocket_sums(idx) = 1.0; %TODO: delete after adapting part 3
            
            % integrate monomers
            gelInfo.species.(species.type).monomer_sums(pos) = 0.0;
            profileData.monomer_sums(idx) = 0.0; %TODO: delete after adapting part 3
            
            % intergrate smear
            gelInfo.species.(species.type).smear_sums(pos) = 0.0;
            profileData.smear_sums(idx) = 0.0; %TODO: delete after adapting part 3
        else
            % integrate aggregates
            pocket_boundaries = [max(1,round(mu-sig)) round(mu+sig)];
            pocket_sums = sum(species.fullprofiles{pos}(pocket_boundaries(1):pocket_boundaries(2) ));
            gelInfo.species.(species.type).pocket_sums(pos) = pocket_sums;
            profileData.pocket_sums(idx) = pocket_sums; %TODO: delete after adapting part 3
            
            
            % integrate monomers
            mu_mono = species.positions(pos);
            sig_mono = species.band_width(pos) * gelInfo.sigma_integrate;
            monomer_boundaries = [max(1,round(mu_mono-sig_mono)) round(mu_mono+sig_mono)];
            monomer_sums = sum(species.fullprofiles{pos}(monomer_boundaries(1):monomer_boundaries(2)));
            gelInfo.species.(species.type).monomer_sums(pos) = monomer_sums;
            profileData.monomer_sums(idx) = monomer_sums; %TODO: delete after adapting part 3
            
            % intergrate smear
            smear_boundaries = [pocket_boundaries(2) monomer_boundaries(1)];
            smear_sums = sum(species.fullprofiles{pos}(smear_boundaries(1):smear_boundaries(2)));
            gelInfo.species.(species.type).smear_sums(pos) = smear_sums;
            profileData.smear_sums(idx) = smear_sums; %TODO: delete after adapting part 3
        end
    end
    
    
    

    
end

