function phi = invert_shadowgraphy(imsource, imtarget)
  refresh_interval = 100;

  % normalisation
  addpath('../lib/minFunc_2012');
  addpath('../lib/minFunc_2012/autoDif');
  imsource = imsource / sum(imsource(:)) * prod(size(imsource));
  imtarget = imtarget / sum(imtarget(:)) * prod(size(imsource));

  % initialise u to be no deflection
  [Nx1, Nx2] = size(imsource);
  [x2, x1] = meshgrid([0:Nx2+1], [0:Nx1+1]); % pad by 1 pixel on each end
  u_init = 0.5 * (x1.^2 + x2.^2);
  u0 = u_init(:);

  % get the interpolant
  xx1 = x1(2:end-1,2:end-1);
  xx2 = x2(2:end-1,2:end-1);
  F = griddedInterpolant(xx1, xx2, imtarget, 'linear', 'nearest'); % mean(imsource(:)));

  % get the grid data point on the source plane
  xx1 = x1(2:end-1,2:end-1);
  xx2 = x2(2:end-1,2:end-1);

  % define the objective function and do L-BFGS optimization
  func_obj = @(u) func_grad(F, imsource, u);

  % search the step size
  [f0, du0] = func_obj(u0);
  step0 = 1e0;
  step = step_search(func_obj, u0, f0, du0, step0);
  u = u0 - step * du0;

  % iteration
  niter = 0;
  f = inf;
  while (f > 1e-3)
    [f, du] = func_obj(u);
    if (f > f0)
      step = step / 2;
      u = u0 - step * du0;
    else
      f0 = f;
      u0 = u;
      du0 = du;
      u = u0 - step * du0;
    end
    niter = niter + 1;
    if mod(niter, refresh_interval) == 0
      step = step_search(func_obj, u0, f0, du0, step);
      fprintf('%6d    %.6e     %.3e\n', niter, f, step);
    end
  end

  phi = u_init - reshape(u, size(u_init));
  phi = phi(2:end-1, 2:end-1);
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
  imratio = imsource ./ imtarget_s;
  dudt = zeros(size(imsource)+2);
  dudt(2:end-1,2:end-1) = -log(imsource ./ (imtarget_s .* det_jac));

  % error and the gradient
  f = max(abs(dudt(:)));
  df = -dudt(:);
end

function step = step_search(func_obj, u0, f0, du0, step)
  % perform the line search with golden ratio algorithm

  % get the bounds
  f = func_obj(u0 - step * du0);
  step0 = step * (f < f0);
  golden = (1 + sqrt(5)) / 2;
  while (f < f0)
    step = step * 2;
    f = func_obj(u0 - step * du0);
  end

  % now the search is between [0, step]
  golden_inv = 1/golden;
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
