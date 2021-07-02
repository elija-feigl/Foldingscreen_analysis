function [data_out, cur_fig] = plot_analysis(gelInfo, gelData)   
    
    %TODO: delete after a better way to plot
    
    %% find best folding using migrate quality metric
    n = length(gelInfo.loop_indices);
    sort_idx = sort(gelInfo.loop_indices);
    
    index_Tscrn = find(startsWith(gelInfo.lanes, 'T', 'IgnoreCase',true));
    index_Mgscrn = find(startsWith(gelInfo.lanes, 'M', 'IgnoreCase',true));
    index_RM = find(startsWith(gelInfo.lanes, 'RM', 'IgnoreCase',true));
    
    index_foldings = sort([index_Tscrn index_Mgscrn index_RM]); 
    
    function index = best_lane(metric, indices)
        if ~isempty(indices)
            [~, i_sort] = sort(metric(indices), 'descend');
            index = indices(i_sort(1));
        end
    end

    
    fraction_monomer = zeros(n, 1);
    fraction_smear = zeros(n, 1);
    fraction_pocket = zeros(n, 1);
    folding_quality_metric = zeros(n, 1);
    folding_quality_metric_migrate = zeros(n, 1);
    rel_mono_migrate = zeros(n, 1);
    rel_mono_spread = zeros(n, 1);
    mono_spread = zeros(n, 1);
    mono_migrate = zeros(n, 1);
    band_fits  = zeros(n, 3);
    
    ladder_migrate_error = gelInfo.species.ladder.migrate_error;
    ladder_spread_error = gelInfo.species.ladder.spread_error;
    
    for i = 1:n
        idx = sort_idx(i);
        [spec,  pos] = get_lanes_from_idx(gelInfo, idx);
        
        fraction_monomer(idx) = gelInfo.species.(spec.type).fraction_monomer(pos);
        fraction_smear(idx) = gelInfo.species.(spec.type).fraction_smear(pos);
        fraction_pocket(idx) = gelInfo.species.(spec.type).fraction_pocket(pos);
        folding_quality_metric(idx) = gelInfo.species.(spec.type).folding_quality_metric(pos);
        folding_quality_metric_migrate(idx) = gelInfo.species.(spec.type).folding_quality_metric_migrate(pos);
        rel_mono_migrate(idx) = gelInfo.species.(spec.type).rel_band_migrate(pos);
        rel_mono_spread(idx) = gelInfo.species.(spec.type).rel_band_spread(pos);
        mono_spread(idx) = gelInfo.species.(spec.type).band_spread(pos);
        mono_migrate(idx) = gelInfo.species.(spec.type).migration_distance(pos);
        band_fits(idx, :) = gelInfo.species.(spec.type).fits(pos, :);
    end
    
    
    index_best = best_lane(folding_quality_metric_migrate, index_foldings);
    if ~isempty(index_Tscrn)
        index_best_Tscrn = best_lane(folding_quality_metric_migrate, index_Tscrn);
        data_out.bestTscrn = gelInfo.lanes{index_best_Tscrn};
    else
        index_best_Tscrn = nan;
    end
    
    if ~isempty(index_Mgscrn)
        index_Mgscrn = index_Mgscrn(index_Mgscrn~=0); % remove zeros if people did not include all Mg samples
        index_best_Mgscrn = best_lane(folding_quality_metric_migrate, index_Mgscrn);
        data_out.bestMgscrn = gelInfo.lanes{index_best_Mgscrn};
    else
        index_best_Mgscrn = nan;
    end
    
    if ~isempty(index_RM)
        index_RM = index_RM(index_RM~=0); % remove zeros if people did not include all RM samples
        index_best_RM = best_lane(folding_quality_metric_migrate, index_RM);
        data_out.bestRM = gelInfo.lanes{index_best_RM};
    else
        index_best_RM = nan;
    end
    
    %% save 
    
    data_out.fractionMonomer = fraction_monomer;
    data_out.fractionSmear = fraction_smear;
    data_out.fractionPocket = fraction_pocket;
    data_out.bestFolding = gelInfo.lanes{index_best};
    data_out.bestFoldingIndex = index_best;
    data_out.bestTscrnIndex = index_best_Tscrn ; 
    data_out.bestMgscrnIndex = index_best_Mgscrn;
    data_out.bestRMIndex = index_best_RM;

    %data_out.M20BetterThanTscrn = M20_better_than_bestT;
    data_out.bandWidthNormalized = mono_spread;
    data_out.migrationDistanceNormalized = mono_migrate;
    data_out.qualityMetric = folding_quality_metric;
    
    %% Plot
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 20 30], 'PaperSize', [20 30]);
    subplot(5,1,1:2)
    % AGE gel with band fits (best highlighted)
    plot_band_fits(gelData, gelInfo)
    shift = gelInfo.pocket.bounding_box(1);
    lanes_bounding_box = gelInfo.lanes_bounding_box;
    
    plot(mean([lanes_bounding_box(index_best,1) lanes_bounding_box(index_best, 2) ]) - shift, ...
            band_fits(index_best,2), 'bo');
        
    plot(mean([lanes_bounding_box(index_best_Tscrn,1) lanes_bounding_box(index_best_Tscrn, 2) ]) - shift , ...
            band_fits(index_best_Tscrn,2), 'b+');
        
    plot(mean([lanes_bounding_box(index_best_Mgscrn,1) lanes_bounding_box(index_best_Mgscrn, 2) ]) - shift , ...
            band_fits(index_best_Mgscrn, 2), 'bx');
        
    title(['Band positions with sigma=' num2str(gelInfo.sigma_integrate)]);

    % fraction per lane plot
    subplot(5,1,3)
    plot(fraction_monomer, '.-'), hold on
    plot(fraction_smear, '.-'), hold on
    plot(fraction_pocket, '.-'), hold on
    ylabel('Fraction')
    set(gca, 'XTick', (1:n), 'XTickLabels', gelInfo.lanes, 'XLim', [1 n])
    legend({'monomer', 'smear', 'pocket'}, 'location', 'best')
    grid on

    % migration distance & spread plot
    err_migrate = ladder_migrate_error * ones(length(rel_mono_migrate),1);
    err_spread = ladder_spread_error * ones(length(rel_mono_spread),1);
    subplot(5,1,4)
    errorbar(rel_mono_spread, err_spread, '.-'), hold on
    errorbar(rel_mono_migrate, err_migrate, '.--'), hold on
    ylabel({'Normalized ', 'rel_spread or rela. migr. distance'})
    set(gca, 'XTick', (1:n), 'XTickLabels', gelInfo.lanes, 'XLim', [1 n], 'YLim', [0 rel_mono_migrate(1)])
    legend({'Band spread', 'Migr. distance'}, 'location', 'best')
    grid on 

    % quality metric plot
    subplot(5,1,5)
    plot(folding_quality_metric, '.-'), hold on
    plot(folding_quality_metric_migrate, '.--'), hold on
    plot(index_best, folding_quality_metric_migrate(index_best), 'o'), hold on
    plot(index_best_Tscrn, folding_quality_metric_migrate(index_best_Tscrn), '+'), hold on
    plot(index_best_Mgscrn, folding_quality_metric_migrate(index_best_Mgscrn), 'x'), hold on
    ylabel({'Quality metric ', 'Monomer fraction/norm_width'})
    set(gca, 'XTick', (1:n), 'XTickLabels', gelInfo.lanes, 'XLim', [1 n])
    legend({'Quality metric', 'Quality metric migrate'}, 'location', 'best')
    grid on

end
    