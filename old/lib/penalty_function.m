%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the penalty function and its gradient from Q. Merigot 2011 of weights.
% Input:
%   * p: coordinate of points from source (numPoints x 2)
%   * lambdap: value of each point (numPoints x 1)
%   * targetDensity: map of the density on target plane (Ny x Nx)
%   * weights: weights of each point for the power cell on target plane from the source (numPoints x 1)
%   * verbose: flag to indicate whether to print message and show image or not (default: 0):
%              0 - no message and no processing images shown
%              1 - shows messages, without images
%              2 - shows messages and images
%   * record_to: record the progress to the specified .mat file, leave '' if not recording (default: '')
% Output:
%   * Phi: value of the penalty function (1 x 1)
%   * delPhi: gradient of the penalty function (numPoints x 1)
% Assumptions:
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1
%   * For all i: 1 <= p(i,2) < size(targetDensity,2)+1 and 1 <= p(i,1) < size(targetDensity,1)+1
%   * sum(targetDensity(:)) == sum(lambdap)

function [Phi, delPhi] = penalty_function(p, lambdap, targetDensity, weights, verbose, record_to)
    if nargin < 5
        verbose = 0;
    end
    if nargin < 6
        record_to = '';
    end

    % get the rectangle bound
    [Ny, Nx] = size(targetDensity);
    pRect = [0 0; 0 Ny; Nx Ny; Nx 0] + 1;

    % obtain the polygons of the power diagram
    [V,C] = power_bounded(p(:,1), p(:,2), weights, pRect);
    % PD = powerDiagramRectBounded(p, weights, pRect);
    % PD = PD{1};

    Ip = zeros([length(C),1]);
    Ap = zeros([length(C),1]);

    % iterate for each point p
    % disp(length(PD));
    parfor (i = [1:length(C)])
        if (length(C{i}) == 0) continue; end;

        % polynomial of the power cell
        xPoly = V(C{i},1);
        yPoly = V(C{i},2);
        [xPoly, yPoly] = poly2cw(xPoly, yPoly);

        % position of the point
        xp = p(i,1);
        yp = p(i,2);

        % get the area and inertia on target density of the power cell
        [A, xcm, ycm, I0] = poly_pixel_area_cm_inertia(xPoly, yPoly, targetDensity);
        Ix0 = I0 - A*(xcm^2+ycm^2) + A*((xp-xcm)^2 + (yp-ycm)^2);

        Ap(i) = A;
        Ip(i) = Ix0;
    end

    if (verbose == 2)
        close all;
        figure;
        imagesc([1:Nx]+0.5, [1:Ny]+0.5, targetDensity); hold on;
        plot(p(:,1), p(:,2), 'b.');

        for (i = [1:length(C)])
            if (length(C{i}) == 0) continue; end;
            % polynomial of the power cell
            xPoly = V(C{i},1);
            yPoly = V(C{i},2);
            patch(xPoly, yPoly, [1 1 1], 'FaceAlpha', 0);
        end
        pause(0.1);
    end

    % Phi = weights' * (lambdap + Ap) - sum(Ip);
    % delPhi = lambdap + Ap;
    Phi = weights' * (lambdap - Ap) + sum(Ip);
    delPhi = lambdap - Ap;

    Phi = -Phi;
    delPhi = -delPhi;
    % disp(norm(delPhi));

    % record the progress
    if (~strcmp(record_to, ''))
        load(record_to, 'all_weights', 'all_penalties');
        all_weights(:,end+1) = weights;
        all_penalties(end+1) = Phi;
        save(record_to, 'all_weights', 'all_penalties', '-append');
    end
end
