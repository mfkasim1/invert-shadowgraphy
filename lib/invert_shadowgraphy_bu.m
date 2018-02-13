function phi = invert_shadowgraphy(imsource, imtarget)
  % parameters
  step = 1e-1;
  tol = 1e-3;
  disp_interval = 1000;

  % normalisation
  imsource = imsource / sum(imsource(:));
  imtarget = imtarget / sum(imtarget(:));

  % initialise u to be no deflection
  [Nx1, Nx2] = size(imsource);
  [x2, x1] = meshgrid([0:Nx2+1], [0:Nx1+1]); % pad by 1 pixel on each end
  u0 = 0.5 * (x1.^2 + x2.^2);
  u = u0;

  % get the grid data point on the source plane
  xx1 = x1(2:end-1,2:end-1);
  xx2 = x2(2:end-1,2:end-1);

  % create the gridded interpolant
  g = griddedInterpolant(xx1, xx2, imtarget, 'linear', 'nearest'); % mean(imsource(:)));

  niter = 0;
  init_time = tic;
  while 1
    % compute the shifted versions of the u
    u_i_j     = u(2:end-1,2:end-1);
    u_ip1_j   = u(3:end,2:end-1);
    u_i_jp1   = u(2:end-1,3:end);
    u_im1_j   = u(1:end-2,2:end-1);
    u_i_jm1   = u(2:end-1,1:end-2);
    u_ip1_jp1 = u(3:end,3:end);
    u_im1_jm1 = u(1:end-2,1:end-2);
    u_im1_jp1 = u(1:end-2,3:end);
    u_ip1_jm1 = u(3:end,1:end-2);

    % get the first and second derivatives of u
    ux1   = (u_ip1_j - u_im1_j) / 2;
    ux2   = (u_i_jp1 - u_i_jm1) / 2;
    ux1x1 = u_ip1_j + u_im1_j - 2*u_i_j;
    ux2x2 = u_i_jp1 + u_i_jm1 - 2*u_i_j;
    ux1x2 = (u_ip1_jp1 + u_im1_jm1 - u_im1_jp1 - u_ip1_jm1) / 4;

    % calculate the determinant of the jacobian
    det_jac = abs(ux1x1 .* ux2x2 - ux1x2 .* ux1x2);

    % the dudt needs to be padded by zeros at the boundaries
    dudt = zeros(size(imtarget)+2);

    % get the target plane intensity if brought to the source plane
    imtarget_s = g(ux1, ux2);
    % imtarget_s = interp2(xx2, xx1, imtarget, ux2, ux1, 'linear', mean(imsource(:)));

    % calculate the update
    dudt(2:end-1,2:end-1) = -log(imsource ./ (imtarget_s .* det_jac));

    % get the new u
    u = u + step * dudt;
    niter = niter + 1;

    if mod(niter, disp_interval) == 0
      fprintf('%6d      %.7e      %.3e\n', niter, max(abs(dudt(:))), toc(init_time));
    end
    if max(abs(dudt(:))) < tol
      break;
    end

  end

  phi = u0 - u;
  phi = phi(2:end-1, 2:end-1);
end
