
function plot_band_fits(gelData, profileData, gelInfo)
%% plotting lanes and pocket fits
    img_L = gelInfo.species.pocket.positions(1);
    img_R = img_L + gelInfo.species.pocket.positions(3);
    img_plot = gelData.images{1}(:,img_L : img_R);
    imagesc(img_plot, [0 3.*std(img_plot(:))]), axis image, colormap gray, hold on
    %n = length(profileData.profiles);

    % plot pocket fits
    mu_p = gelInfo.species.pocket.fits(2);
    sig_p = gelInfo.species.pocket.fits(3) * profileData.sigma_integrate;
    %limit = profileData.has_ladder +1;
    
    function get_fits_new(gelInfo, species, profileData, mu_p, img_L, sig_p)
        index = species.index;
        for i=1:length(index)
            x = profileData.lanePositions(index(i),1:2) - double(img_L);
            mu = species.fits(i,2);
            if species.type == "mono"
                st = gelInfo.species.staple.fits(i,2);
            else
                st = [0 0 0];
            end
            sig = species.fits(i,3) * profileData.sigma_integrate;
            
            if species.type == "staple"
                c = [1.0, 0.0, 0.0];  
            
            elseif mu < mu_p
                c = [1.0 1.0 1.0];
            
            else        
                c = [1.0 0.77 0.13];
            end
        plot(mean(x), mu_p, '.', 'color',  c);
        plot([mean(x) mean(x)], [mu_p-sig_p mu_p+sig_p], 'color',  c);
        % plot leading band fits
        plot([x(1) x(2)] , [mu mu], 'color',  c);
        plot([mean(x) mean(x)], [mu-sig mu+sig],  'color', c);
        plot(mean(x), st, '.', 'color', [1.0 0.77 0.13]);
        end
    end

    if gelInfo.species.ladder.type ~= "none"
        get_fits_new(gelInfo, gelInfo.species.ladder, profileData, mu_p, img_L, sig_p)
    end
    
    if gelInfo.species.scaffold.type ~= "none"
        get_fits_new(gelInfo, gelInfo.species.scaffold, profileData, mu_p, img_L, sig_p)
    end
    get_fits_new(gelInfo, gelInfo.species.mono, profileData, mu_p, img_L, sig_p)
        

end

