verbose = 1;
scale = 4;
N = 50*scale; % image size: N x N pixels (default: 50)
sourceMap = ones(N); % uniform source distribution profile

% construct the deflection potential
[X,Y] = meshgrid(linspace(-1, 1, N)); % create a coordinate with N points going from -1 to 1
sigma = 0.25; % sigma of the gaussian profile
Phi0 = 60*scale*scale; % peak value of the deflection potential (default: 80)
Phi = Phi0 * exp(-(X.^2 + Y.^2)/2/sigma^2);

% get the shadowgraph image with the corresponding deflection potential
if (verbose) disp('Obtaining the shadowgram ...'); end;
forwardTic = tic;
targetMap = forward_shadowgraphy(sourceMap, Phi); % usually it takes ~2s in my computer
forwardTime = toc(forwardTic);
if (verbose) disp(sprintf('Finish in %fs', forwardTime)); end

imagesc(targetMap);
