function [imtarget] = forward_shadowgraphy(imsource, phi)
  % calculate the total deflection, so that y = del(u)
  [Nx1, Nx2] = size(imsource);
  [x2, x1] = meshgrid([0:Nx2+1], [0:Nx1+1]); % pad by 1 pixel on each end
  phi_pad = [phi(1,:);phi;phi(end,:)];
  phi_pad = [phi_pad(:,1),phi_pad,phi_pad(:,end)];
  u = 0.5 * (x1.^2 + x2.^2) - phi_pad;

  % compute the shifted versions of the u
  u_i_j = u(2:end-1,2:end-1);
  u_ip1_j = u(3:end,2:end-1);
  u_i_jp1 = u(2:end-1,3:end);
  u_im1_j = u(1:end-2,2:end-1);
  u_i_jm1 = u(2:end-1,1:end-2);
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

  % calculate the map from the source to the target plane
  y1 = ux1;
  y2 = ux2;

  % interpolate the determinant and the source plane
  xx1 = x1(2:end-1,2:end-1);
  xx2 = x2(2:end-1,2:end-1);
  imsource_t = griddata(y2, y1, imsource, xx2, xx1);
  det_jac_t  = griddata(y2, y1, det_jac , xx2, xx1);
  imtarget   = imsource_t ./ det_jac_t;
end
