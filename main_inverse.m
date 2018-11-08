function phi = main_inverse(imsource, imtarget, options)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Retrieving the deflection potential given the source image (i.e. the image
  % without the deflection) and the target image (i.e. the image with the
  % deflection) using the algorithm suggested in Sulman, et al. (2011), "An
  % efficient approach for the numerical solution of the Monge–Ampère equation"
  % combined with adaptive step search.
  %
  % Input:
  % * imsource: the source image, for undefined region, can be set to nan
  % * imtarget: the target image, for undefined region, can be set to nan
  %
  % Output:
  % * phi: the deflection potential
  %
  % Constraints:
  % * imsource and imtarget need no normalization, it normalizes inside the
  %     function to have mean = 1
  % * the defined region must have the same shape and position
  % * the values uses the normalization values where each pixel has size of 1
  % * the defined region should be convex
  %
  % Author: Muhammad F. Kasim (University of Oxford, 2018)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (nargin < 3)
      options.interp = 'linear';
  end

  % options
  defopt.refresh_interval = 100;
  defopt.rel_tol = 1e-3;
  defopt.max_niter = 50000;
  defopt.minstep = 1e-3;
  defopt.max_time = 60*5;
  defopt.max_niter_no_update = 5000;
  defopt.alpha = 0.1;
  % it seems that nearest neighbor interpolation is more robust than 'linear'
  defopt.interp = 'nearest';
  defopt.extrap = 'nearest';

  % get the options and default options field names
  defopt_fields = fieldnames(defopt);
  opt_fields = fieldnames(options);

  % if there is a defopt field that is not in opt field, set opt field to default
  for idef = 1:length(defopt_fields)
    defopt_field = defopt_fields{idef};
    found = false;

    for io = 1:length(opt_fields)
      opt_field = opt_fields{io};
      if strcmp(opt_field, defopt_field)
        found = true;
        break;
      end
    end

    % if not found, then set the default value
    if (~found)
      options = setfield(options, defopt_field, getfield(defopt, defopt_field));
    end
  end

  % set the parameter names
  refresh_interval = options.refresh_interval;
  rel_tol = options.rel_tol;
  max_niter = options.max_niter;
  minstep = options.minstep;
  max_time = options.max_time;
  max_niter_no_update = options.max_niter_no_update;
  alpha = options.alpha;
  interp = options.interp;
  extrap = options.extrap;

  % get the mask
  mask = ~isnan(imtarget);

  % normalisation
  imsource = imsource / mean(imsource(mask));
  imtarget = imtarget / mean(imtarget(mask));

  % extrapolate outside the mask
  imtarget_fill = fillmissing(imtarget, 'constant', nan);

  % initialise u to be no deflection
  [Nx1, Nx2] = size(imsource);
  [x2, x1] = meshgrid([0:Nx2+1], [0:Nx1+1]); % pad by 1 pixel on each end
  u_init = 0.5 * (x1.*x1 + x2.*x2);
  u0 = u_init(:);

  % get the interpolant
  xx1 = x1(2:end-1,2:end-1);
  xx2 = x2(2:end-1,2:end-1);
  F = griddedInterpolant(xx1, xx2, imtarget_fill, interp, extrap);

  % get the grid data point on the source plane
  xx1 = x1(2:end-1,2:end-1);
  xx2 = x2(2:end-1,2:end-1);

  % define the objective function and do L-BFGS optimization
  func_obj = @(u) func_grad(F, imsource, u);

  % minimize the objective
  options.max_niter = max_niter;
  options.rel_tol = rel_tol;
  options.max_time = max_time;
  options.max_niter_no_update = max_niter_no_update;
  options.minstep = minstep;
  options.refresh_interval = refresh_interval;
  options.alpha = alpha;
  [umin, fmin] = minimize(func_obj, u0, options);
  % opts.MaxFunEvals = 50000;
  % opts.Method = 'csd';
  % [umin, fmin] = minFunc(func_obj, u0, opts);

  phi = u_init - reshape(umin, size(u_init));
  phi = phi(2:end-1, 2:end-1);
end

function [umin, fmin] = minimize(func_obj, u0, options)
  minstep = options.minstep;
  alpha = options.alpha;
  stime = tic;

  % search the step size
  [f0, du0] = func_obj(u0);
  finit = f0;
  step0 = minstep;
  step = step_search(func_obj, u0, f0, du0, step0);
  u = u0 - step * du0;

  % iteration
  niter = 0;
  f = inf;
  fmin = inf;
  umin = 0;
  last_niter_update = 0;
  while 1
    [f, du] = func_obj(u);
    if f < fmin
      fmin = f;
      umin = u;
      last_niter_update = niter;
    end
    if f > f0 && step > minstep
      step = step / 2;
      u = u0 - step * du0;
    else
      f0 = f;
      u0 = u;
      du0 = du * alpha + du0 * (1-alpha);
      u = u0 - step * du0;
    end
    niter = niter + 1;
    if step < minstep
      step = minstep;
    end
    if mod(niter, options.refresh_interval) == 0 || niter == 1
      step = step_search(func_obj, u0, f0, du0, step);
      if step < minstep
        step = minstep;
      end
      fprintf('%6d    %.6e     %.3e    %.3es\n', niter, f, step, toc(stime));
    end

    % stopping conditions
    if niter > options.max_niter
      break;
    end
    if fmin / finit < options.rel_tol
      break;
    end
    if toc(stime) > options.max_time
      break;
    end
    if niter - last_niter_update > options.max_niter_no_update
      break
    end
  end
  fprintf('%6d    %.6e     %.3e    %.3es\n', niter, f, step, toc(stime));
end

function [f,df] = func_grad(F, imsource, u)
  % resize u
  [Nx1, Nx2] = size(imsource);
  u = reshape(u, (Nx1+2), (Nx2+2));

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
  ux1 = (u_ip1_j - u_im1_j) / 2;
  ux2 = (u_i_jp1 - u_i_jm1) / 2;
  ux1x1 = u_ip1_j + u_im1_j - 2*u_i_j;
  ux2x2 = u_i_jp1 + u_i_jm1 - 2*u_i_j;
  ux1x2 = (u_ip1_jp1 + u_im1_jm1 - u_im1_jp1 - u_ip1_jm1) / 4;

  % calculate the determinant of the jacobian
  det_jac = abs(ux1x1 .* ux2x2 - ux1x2 .* ux1x2);

  % get the target plane intensity if brought to the source plane
  imtarget_s = F(ux1, ux2);
  dudt = zeros(size(imsource)+2);
  dudt(2:end-1,2:end-1) = -log(abs(imsource ./ (imtarget_s .* det_jac)));

  % normalise the not-normal values
  dudt(isnan(dudt)) = 0;
  dudt(isinf(dudt)) = 0; %max(dudt(~isinf(dudt)));

  % error and the gradient
  f = mean((dudt(:)).^2);
  df = -dudt(:);
end

function step = step_search(func_obj, u0, f0, du0, step)
  % perform the line search with golden ratio algorithm
  golden = (1 + sqrt(5)) / 2;
  golden_inv = 1/golden;

  % get the bounds
  f = func_obj(u0 - step * du0);
  step0 = step * (f < f0);
  while (f < f0)
    step = step * 2;
    f = func_obj(u0 - step * du0);
  end

  % now the search is between [0, step]
  a = step0;
  b = step;
  c = b - (b - a) * golden_inv;
  d = a + (b - a) * golden_inv;
  fa = f0;
  fb = f;
  fc = func_obj(u0 - c * du0);
  fd = func_obj(u0 - d * du0);
  for i = [1:4]
    if fc < fd
      b = d;
      fb = fd;
    else
      a = c;
      fa = fc;
    end
    c = b - (b - a) * golden_inv;
    d = a + (b - a) * golden_inv;
    if fc < fd
      fd = fc;
      fc = func_obj(u0 - c * du0);
    else
      fc = fd;
      fd = func_obj(u0 - d * du0);
    end
  end
  if fc < fd
    step = c;
  else
    step = d;
  end
end
