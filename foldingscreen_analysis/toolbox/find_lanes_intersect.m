function [ lanes ] = find_lanes_intersect( image, pos )
%% find lanes edges at intersection with user-supplied threshold
%  This function determines the positions of lanes based on a user defined
%  value.
%  Input: image = Matrix of gray values, pos = positions of outer rectangle as [x0 y0 width height]
%  x0, y0 are the top left positions
%  Output: lanes = positions of lanes
%  Example: find_lanes_intersect(some_image, [100 150 200 50 ])

    area = image( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)); % get the image inside pos
    x = double(pos(1):pos(1)+pos(3)); % define x-axis (horizontal axis)

    horizontalProfile = sum(area); % integrate along vertical (y-axis) 

    % plot subimage
    cur_fig = figure;


    %plot horizontal profile and determine threshold interactively
    plot(x, horizontalProfile, 'b');
    hold on
    hline(min(horizontalProfile), 'k--');
    hline(max(horizontalProfile), 'g--');
    threshold_init = (max(horizontalProfile)-min(horizontalProfile))/4; % initialize threshold
    hline(threshold_init, 'r--');
    set(gca, 'XLim', [pos(2) pos(2)+pos(4)])
    plot(x, horizontalProfile, 'b')
    xlim = [x(1) x(end)]; 
    set(gca, 'XLim', xlim)
    xlabel('Horizontal Position [pixel]')
    ylabel('Intensity [a.u.]')
    
    % create draggable line and read the choosen value
    h = imline(gca, xlim, [threshold_init threshold_init]); %plot horzontal line
    setColor(h,[1 0 0]);
    setPositionConstraintFcn(h, @(pos)[  xlim' [pos(1,2);pos(1,2)]  ])
    pos_line = wait(h);
    threshold = pos_line(1,2);

    close(cur_fig)


    % deterime intersection of lines by first subtracting the threshold and
    % then finding roots
    y_shift = horizontalProfile - threshold;

    %find roots of y_shift based on change of sign
    r = [];
    for i=1:length(y_shift)-1
        if y_shift(i)*y_shift(i+1) < 0 %root
           r = [r ; i];
        end
    end
    
    
    % check if an even number of roots where found
    if mod(length(r),2) == 1 % uneven 
        cur_fig = figure();

        % plot found lanes
        subplot(2, 1, 1)
        imagesc(image), axis image, colormap gray
        set(gca, 'XLim', [pos(1) pos(1)+pos(3)], 'YLim', [pos(2) pos(2)+pos(4)])
        
        % plot profile
        subplot(2, 1, 2)
        h_leg(1) = plot(x, horizontalProfile, 'b');
        h_leg(2:1+length(r)) = vline(x(r), 'r');
        legend(h_leg(1:2), {'Horizontal Profile', 'Found Intersections'})
        set(gca, 'XLim', [pos(1) pos(1)+pos(3)])
        for i=1:length(r)
            text(x(r(i)), threshold, num2str(i) ) % plot index of intersection
        end

        answer = inputdlg({'Uneven number of intersections. Which intersection should be dismissed:'},'Dismiss intersect',1,{'1'});
        dismiss_i = str2double(answer{1});
        r(dismiss_i) = [];
        
        close(cur_fig)
    end
    
    % write output
    lanes = zeros(length(r)/2,4);
    for i=1:size(lanes, 1)
       lanes(i, 2)= pos(2); % top y-positions stays constant
       lanes(i, 4)= pos(4); % height stays constant
       lanes(i, 1) = x(r(i*2-1)); % top x-position
       lanes(i, 3) = x(r(i*2))-x(r(i*2-1)); % width
    end

end

