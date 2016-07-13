%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to visualise pixels, lines, and polygon
% Input:
%   * x0, y0: x and y-coordinates of the points (1xN)
%   * drawType: can be 'pixels', 'lines', or 'polygon'
%   * color: the first letter of the color in string

function visualise_pixels_polygon(x0, y0, drawType, color)
    if (strcmp(drawType, 'pixels'))
        hold on;
        for (i = [1:length(x0)])
            xp = [0,0,1,1,0] + x0(i);
            yp = [0,1,1,0,0] + y0(i);
            
            patch(xp, yp, 1, 'FaceColor', color);
        end
    elseif (strcmp(drawType, 'lines'))
        hold on;
        for (i = [1:length(x0)-1])
            plot(x0(i:i+1), y0(i:i+1), color);
        end
    elseif (strcmp(drawType, 'polygon'))
        hold on;
        x0(end+1) = x0(1);
        y0(end+1) = y0(1);
        for (i = [1:length(x0)-1])
            plot(x0(i:i+1), y0(i:i+1), color);
        end
    end
end
