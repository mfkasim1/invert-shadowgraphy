%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predict the intensity on the target plane given the deflection potential, Phi, and
%   the intensity profile on the source plane, sourceMap.
% The beam from source plane is mapped to the target plane as function of its position on the source plane.
% The map is given by:
%   xt(xs,ys) = xs - \partial(Phi(xs,ys))/\partial(xs);
%   yt(xs,ys) = ys - \partial(Phi(xs,ys))/\partial(ys);
% Input:
%   * sourceMap: intensity profile on the source plane (Ny x Nx)
%   * Phi: deflection potential according to the equation above (Ny x Nx)
% Output:
%   * targetMap: intensity profile on the target plane (Ny x Nx)
% Assumptions:
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1
%   * sourceMap and Phi must have the same size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targetMap = main_forward(sourceMap, Phi)
    add_libs;

    % getting the basic parameters for this algorithm
    [Ny, Nx] = size(sourceMap);
    [x,y] = meshgrid([1:Nx], [1:Ny]);
    maxIntensity = max(sourceMap(:)) * 1000;

    % obtain the deflection angles
    alphax = -diff(Phi, [], 2);
    alphay = -diff(Phi, [], 1);
    alphax(:,Nx) = alphax(:,Nx-1); % just to make it the same size
    alphay(Ny,:) = alphay(Ny-1,:);

    % get the mapping plan
    xmap = x + alphax;
    ymap = y + alphay;

    I = zeros(size(x));
    times = zeros([1 7]); % optim:

    % start the algorithm here
    tic;
    [Nymap, Nxmap] = size(xmap);
    xA = xmap(1:Nymap-1, 1:Nxmap-1); yA = ymap(1:Nymap-1, 1:Nxmap-1);
    xB = xmap(2:Nymap  , 1:Nxmap-1); yB = ymap(2:Nymap  , 1:Nxmap-1);
    xC = xmap(2:Nymap  , 2:Nxmap  ); yC = ymap(2:Nymap  , 2:Nxmap  );
    xD = xmap(1:Nymap-1, 2:Nxmap  ); yD = ymap(1:Nymap-1, 2:Nxmap  );

    % sort polygons to have cw representations
    Ntot = prod(size(xA));
    xpoly = [reshape(xA, [1, Ntot]);     ...
             reshape(xB, [1, Ntot]); ...
             reshape(xC, [1, Ntot]); ...
             reshape(xD, [1, Ntot])];
    ypoly = [reshape(yA, [1, Ntot]); ...
             reshape(yB, [1, Ntot]); ...
             reshape(yC, [1, Ntot]); ...
             reshape(yD, [1, Ntot])];

    xcentre = repmat(mean(xpoly,1), [4,1]);
    ycentre = repmat(mean(ypoly,1), [4,1]);
    dxpoly = xpoly - xcentre;
    dypoly = ypoly - ycentre;
    anglepoly = atan2(dypoly, dxpoly);
    [~, idxsort] = sort(anglepoly, 1, 'descend');
    idxsort = idxsort + size(idxsort,1) * repmat([0:size(idxsort,2)-1], [size(idxsort,1), 1]);
    xparallelogram = xpoly(idxsort);
    yparallelogram = ypoly(idxsort);

    % get the area of the parallelograms in format:
    % [dS(0); dS(1); ...]
    dS = polyarea(xparallelogram, yparallelogram)';

    % get the normalised coordinate in format:
    % [xA0, xB0, xC0, xD0; xA1, xB1, xC1, xD1; ...];
    ixparallelogram = (xparallelogram - x(1,1))' + 1;
    iyparallelogram = (yparallelogram - y(1,1))' + 1;

    % get the I0 in form of [I0linear(1); I0linear(2), ...];
    I0linear = reshape(sourceMap(1:Nymap-1, 1:Nxmap-1), [(Nymap-1)*(Nxmap-1), 1]);
    I0centre = I0linear ./ dS;
    dt = toc;
    times(1) = times(1) + dt;

    for (i = [1:(Nymap-1)*(Nxmap-1)])
        % if (mod(i-1, (Nxmap-1)) == 0) disp((i-1) / (Nxmap-1) + 2); end;

        tic;
        % get the pixels enclosed by the polygon
        [ixl, iyl, ixp, iyp] = pixels_enclosed_by_polygon(ixparallelogram(i,:), iyparallelogram(i,:));
        dt = toc;
        times(2) = times(2) + dt;

        % iterate for all pixels enclosed by polygon
        tic;
        for (ip = [1:length(ixp)])
            ix = ixp(ip);
            iy = iyp(ip);

            % check the range
            if ((iy < 1) || (iy > size(y,1))) continue; end;
            if ((ix < 1) || (ix > size(x,2))) continue; end;

            I(iy,ix) = I(iy,ix) + I0centre(i);

        end
        dt = toc;
        times(3) = times(3) + dt;

        % iterate for all pixels crossed by polygon's edges
        for (il = [1:length(ixl)])
            tic;
            ix = ixl(il);
            iy = iyl(il);

            % check the range
            if ((iy < 1) || (iy > size(y,1))) continue; end;
            if ((ix < 1) || (ix > size(x,2))) continue; end;
            % construct the voxel x and y-coordinate
            xvox = x(1,ix) + [0;0;1;1];
            yvox = y(iy,1) + [0;1;1;0];
            dt = toc;
            times(4) = times(4) + dt;

            if (dS ~= 0)
                tic;
                [xpoly, ypoly] = clip_polygons_with_rect(xvox, yvox, xparallelogram(:,i), yparallelogram(:,i));
                dt = toc;
                times(5) = times(5) + dt;
                tic;
                % disp(size(xpoly));
                dSIntersect = polyareaconvex(xpoly, ypoly);
                % disp(size(dSIntersect));
                I(iy,ix) = I(iy,ix) + I0centre(i) * dSIntersect;
                dt = toc;
                times(6) = times(6) + dt;
            else
                I(iy,ix) = maxIntensity;
            end
        end
    end

    % clip the intensity
    I(I > maxIntensity) = maxIntensity;

    % copy the intensity at the boundary
    I(:,Nx) = I(:,Nx-1);
    I(Ny,:) = I(Ny-1,:);

    % assign the output
    targetMap = I;

    % I = times; % optim
end
