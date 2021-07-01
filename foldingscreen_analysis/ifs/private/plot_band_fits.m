
function plot_band_fits(gelData, gelInfo)
% plotting lanes and pocket fits

    img_L = gelInfo.pocket.bounding_box(1);
    img_R = img_L + gelInfo.pocket.bounding_box(3);
    img_plot = gelData.images{1}(:,img_L : img_R);
    imagesc(img_plot, [0 3.*std(img_plot(:))]), axis image, colormap gray, hold on


    % plot pocket fits
    mu_p = gelInfo.pocket.position;
    sig_p = gelInfo.pocket.width * gelInfo.sigma_integrate;

    
    function plot_fits(gelInfo, species, pos, idx, mu_p, img_L, sig_p)

        x = gelInfo.lane_bounding_box(idx,1:2) - double(img_L);
        mu = species.band_positions(pos);
        if species.type == "mono"
            st = species.staple.fits(pos,2);
        else
            st = [0 0 0];
        end
        sig = species.band_width(pos) * gelInfo.sigma_integrate;

        if species.type == "mono"
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
    
    for i = 1:length(gelInfo.loop_indices)
        idx = gelInfo.loop_indices(i);

        [species,  pos] = get_lanes_from_idx(gelInfo, idx);
        plot_fits(gelInfo, species, pos, idx, mu_p, img_L, sig_p)
    end   
end

