% NOTE: modified from MATLABTOOLBOX

function [ lanes ] = manual_lane_selection( image, pos, varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
area = image( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)); % get the image inside pos
    x_axis = double(pos(1):pos(1)+pos(3)); % define x-axis (horizontal axis)

    horizontalProfile = sum(area); % integrate along vertical (y-axis) 


    if isempty(varargin)     
        tmp = inputdlg('Number of lanes?', 'Number of lanes', 1, {'1'});
        N_lanes = str2double(tmp{1});
    else
        if varargin{1}>0
            N_lanes = varargin{1};
        else
            tmp = inputdlg('Number of lanes?', 'Number of lanes', 1, {'1'});
            N_lanes = str2double(tmp{1});
        end
    end
            
    cf = figure;
    hold all
    plot(x_axis, horizontalProfile)
    set(gca, 'Xlim', [min(x_axis), max(x_axis)])
    
    x = zeros(N_lanes);
    title('Select intial lanes.')
    for i=1:length(x)
        [x(i), ~] = ginput(1);
        vline(x(i), 'r');
    end
    
    close(cf)

    % write areas
    lanes = zeros(N_lanes,4);
    width = (x(2)-x(1)) * 0.6; 
    for i=1:size(lanes, 1)
       lanes(i, 2)= pos(2); % top y-positions stays constant
       lanes(i, 4)= pos(4); % height stays constant

       lanes(i, 3) = width;
       left_edge = x(i) - 0.5 * width;
       left_edge = max(min(x_axis), left_edge);
       left_edge = min(max(x_axis), left_edge);
       lanes(i, 1) = left_edge;w
    end
   
end

