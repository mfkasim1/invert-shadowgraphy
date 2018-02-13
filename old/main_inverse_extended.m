%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the deflection potential from the sourceMap intensity, given the targetMap intensity.
% The beam from source plane is mapped to the target plane as function of its position on the source plane.
% The map is given by:
%   xt(xs,ys) = xs - \partial(Phi(xs,ys))/\partial(xs);
%   yt(xs,ys) = ys - \partial(Phi(xs,ys))/\partial(ys);
% Input:
%   * sourceMap: intensity profile on the source plane (Ny x Nx)
%   * targetMap: intensity profile on the target plane (Ny x Nx)
%   * num_sites: number of sites (default: floor(0.8 * numel(targetMap)))
%   * algorithm: 'lbfgs' or 'quasi-newton' (default: 'lbfgs' if numel <= 4x10^4, else 'quasi-newton')
%   * verbose: flag to indicate whether to print message and show image or not (default: 0):
%              0 - no message and no processing images shown
%              1 - shows messages, without images
%              2 - shows messages and images
%   * record_to: record the progress to the specified .mat file, leave '' if not recording (default: '')
% Output:
%   * Phi: deflection potential according to the equation above (Ny x Nx)
%   * sites: sampling sites from the sourceMap using Lloyd's algorithm (numPoints x 3)
%   * w: weights for each site in the resulting power diagram (numPoints x 1)
% Assumptions:
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1
%   * sourcemap and targetMap must have the same size
%   * total intensity on the sourceMap and targetMap must equal to 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Phi, sites, w] = main_inverse_extended(sourceMap, targetMap, num_sites, algorithm, verbose, record_to)
    add_libs;
    
    % default parameters
    if nargin < 3
        num_sites = floor(0.8 * numel(targetMap));
    end
    if nargin < 4
        if (numel(targetMap) <= 40000)
            algorithm = 'quasi-newton';
        else
            algorithm = 'lbfgs';
        end;
    end
    if nargin < 5
        verbose = 0;
    end
    if nargin < 6
        record_to = '';
    end

    % getting the basic parameters for this algorithm
    [Ny, Nx] = size(targetMap);
    [X,Y] = meshgrid([1:Nx], [1:Ny]);
    N = num_sites;
    pRect = [0 0; 0 Ny; Nx Ny; Nx 0] + 1;

    % get the sites on the sourceMap
    % sample the source map using rejection algorithm and Lloyd's algorithm
    if (verbose) disp('Sampling the source map'); end
    lloydsOptions = NaN;
    [Px0, Py0] = initial_random_sample(N, sourceMap);
    [Px , Py , Ap] = weighted_lloyds_algorithm(Px0, Py0, sourceMap, lloydsOptions);
    p = [Px, Py];
    lambdap = Ap;

    % save it to sites so people can reuse it for other targetMap (with the same sourceMap)
    sites = zeros([size(p,1),3]);
    sites(:,1:2) = p;
    sites(:,3) = lambdap;

    if (~strcmp(record_to, ''))
        all_weights = [];
        all_penalties = [];
        save(record_to, 'sites', 'sourceMap', 'targetMap', 'all_weights', 'all_penalties');
    end

    % use gradient descent to determine the optimal weight for optimal transport map
    if (verbose) disp('Getting the optimal transport map'); end;
    fun = @(w) penalty_function(p, lambdap, targetMap, w, verbose, record_to);
    w0 = zeros([N,1]);
    max_iter = max(Nx, Ny) * 4;
    if (strcmp(algorithm, 'quasi-newton'))
        if (verbose) displayOpt = 'iter';
        else displayOpt = 'none'; end
        options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'GradObj', 'on', 'Display', displayOpt, 'MaxIter', max_iter); % quasi-newton
        w = fminunc(fun, w0, options);
    elseif (strcmp(algorithm, 'lbfgs'))
        if (verbose) options.Display = 1;
        else options.Display = 0; end;
        options.MaxIter = max_iter;
        w = minFunc(fun, w0, options); % l-bfgs
    end

    % determine the weighted centroid of each power cell
    [V,C] = power_bounded(p(:,1), p(:,2), w, pRect);
    % PD = powerDiagramRectBounded(p, w, pRect);
    % PD = PD{1};
    pt = zeros([length(C),2]);
    numEmpty = 0;
    emptySites = [];
    for (j = [1:length(C)])
        % logging the error
        if (length(C{j}) == 0)
            numEmpty = numEmpty + 1;
            emptySites(numEmpty,:) = p(j,:);
            continue;
        end
        xPoly = V(C{j},1);
        yPoly = V(C{j},2);
        [xPoly, yPoly] = poly2cw(xPoly, yPoly);
        [~, pt(j,1), pt(j,2), ~] = poly_pixel_area_cm_inertia(xPoly, yPoly, targetMap);
    end

    % debugging
    if (numEmpty > 0)
        disp(sprintf('Warning: there are %d sites without cells, and they are at:', numEmpty));
        disp(emptySites);
    end


    % move the points closest to the corners to the corner (to make the convex hull covers all the area)
    ps = p;
    for (j = [1:size(pRect,1)])
        corner = pRect(j,:);
        % determine the point closest to the corner
        vec = bsxfun(@minus, p, corner);
        distance = sum(vec.^2, 2);
        minIdx = find(distance == min(distance), 1);

        % move the point to the corner
        ps(minIdx,:) = corner;
    end

    % interpolate and extrapolate for xMap and yMap
    if (verbose) disp('Interpolating the transportation map'); end
    Fx = scatteredInterpolant(ps(:,1), ps(:,2), pt(:,1), 'natural', 'linear');
    Fy = scatteredInterpolant(ps(:,1), ps(:,2), pt(:,2), 'natural', 'linear');
    xMap = Fx(X+0.5,Y+0.5);
    yMap = Fy(X+0.5,Y+0.5);

    % get the deflection angle
    if (verbose) disp('Calculating the deflection potential'); end
    alphax = xMap - (X+0.5);
    alphay = yMap - (Y+0.5);

    % invert the grad by choosing one point as zero (pZero)
    pZero = floor(size(targetMap)/2);
    cumsumX = cumsum(alphax, 2);
    cumsumY = cumsum(alphay, 1);
    cumsumX = bsxfun(@minus, cumsumX, cumsumX(:,pZero(2)));
    cumsumY = bsxfun(@minus, cumsumY, cumsumY(pZero(1),:));
    PhiX = -bsxfun(@plus, cumsumX, cumsumY(:,pZero(2)));
    PhiY = -bsxfun(@plus, cumsumY, cumsumX(pZero(1),:));
    Phi = (PhiX + PhiY)/2;
    if (verbose) disp('Done!'); end;
end
