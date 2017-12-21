% This file is to demonstrate the forward and invert shadowgraphy
% Author: Muhammad Firmansyah Kasim (University of Oxford, 2016)


% test = 0 to take shadowgraph image of a Gaussian deflection potential (takes < 1 minute in my computer)
% test = 1 to take shadowgraph image of a Gaussian deflection potential and invert the image to retrieve the deflection potential (takes ~10 minutes in my computer)
% test = 2 to demonstrate the use of invert_shadowgraphy function (takes ~10 minutes in my computer)
test = 1;
add_libs;

if test <= 1
    verbose = 1;
    N = 50; % image size: N x N pixels (default: 50)
    sourceMap = ones(N); % uniform source distribution profile

    % construct the deflection potential
    [X,Y] = meshgrid(linspace(-1, 1, N)); % create a coordinate with N points going from -1 to 1
    sigma = 0.25; % sigma of the gaussian profile
    Phi0 = 80; % peak value of the deflection potential (default: 80)
    Phi = Phi0 * exp(-(X.^2 + Y.^2)/2/sigma^2);

    % get the shadowgraph image with the corresponding deflection potential
    if (verbose) disp('Obtaining the shadowgram ...'); end;
    forwardTic = tic;
    targetMap = main_forward(sourceMap, Phi); % usually it takes ~2s in my computer
    forwardTime = toc(forwardTic);
    if (verbose) disp(sprintf('Finish in %fs', forwardTime)); end

    % displaying the source map, deflection potential, and the shadowgraphy image
    close all;
    subplot(2,2,1);
    imagesc(Phi); colormap default; colorbar;
    title('Deflection potential');
    subplot(2,2,2);
    imagesc(targetMap); colormap gray;
    title('Shadowgram image');
    subplot(2,2,4);
    plot(targetMap(ceil(end/2),:));
    title('Central horizontal slice');
end
if test == 1
    % set the inversion algorithm parameters
    num_sites = floor(numel(targetMap) * 0.8); % number of sites to be tried (I recommend 0.8 * number of pixels)
    algorithm = 'quasi-newton'; % 'lbfgs' uses much less memory, 'quasi-newton' gives slightly better performance (for small input size, use quasi-newton)

    % now invert the image
    if (verbose) disp('Retrieving the deflection potential ...'); end
    inverseTic = tic;
    [PhiI, sites, w] = main_inverse_extended(sourceMap, targetMap, num_sites, algorithm, verbose); % it takes ~7 minutes in my computer
    inverseTime = toc(inverseTic);
    if (verbose) disp(sprintf('Finish in %fs', inverseTime)); end

    % normalise the retrieved potential
    PhiI = PhiI - min(PhiI(:));

    % displaying the retrieved deflection potential
    subplot(2,2,1);
    imagesc(PhiI); colormap default; colorbar;
    title('Retrieved deflection potential');
    subplot(2,2,3);
    plot(Phi(ceil(end/2),:), 'b-'); hold on;
    plot(PhiI(ceil(end/2),:), 'g--'); hold off;
    title('Horizontal slice of the deflection potentials');
    subplot(2,2,2);
    imagesc(targetMap); colormap gray;
    title('Shadowgram image');
    subplot(2,2,4);
    plot(targetMap(ceil(end/2),:));
    title('Central horizontal slice');

elseif test == 2
    % get the current directory
    currentDir = mfilename('fullpath');

    % normalise the slash
    currentDir = strrep(currentDir, '\', '/');

    % parse the directory only
    kkkkk = strfind(currentDir, '/');
    currentDir = currentDir(1:kkkkk(end));
    if ~strcmp(currentDir(end), '/') currentDir = strcat(currentDir, '/'); end

    % get the file directory
    fdir = strcat(currentDir, 'test-figures/');
    filename = strcat(fdir, 'test.png');

    % invert it using invert_shadowgraphy
    invert_shadowgraphy(filename, 'num_sites', 2000, 'fdir_out', fdir, 'algorithm', 'lbfgs', 'num_workers', 4);

    % load the result file
    load(strcat(fdir, 'test-data.mat'));

    % display the results
    PhiI = PhiI - min(PhiI(:));
    subplot(2,2,1);
    imagesc(PhiI); colormap default; colorbar;
    title('Retrieved deflection potential');
    subplot(2,2,2);
    imagesc(targetMap); colormap gray;
    title('Shadowgram image');
    subplot(2,2,3);
    plot(PhiI(ceil(end/2),:));
    title('Horizontal slice of the potential');
    subplot(2,2,4);
    imagesc(targetMapI); colormap gray;
    title('Shadowgram image from the retrieved potential');

end
