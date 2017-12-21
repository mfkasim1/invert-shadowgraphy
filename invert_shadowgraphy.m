% Invert the grayscale shadowgraphy image
% Input:
%   * filename: string to a filename
% Options:
%   * 'num_sites': number of sites (default: min(100000, 0.8*number of pixels of the image))
%   * 'source_map': 0 - uniform with the same size as the file (default)
%                   >0 - using tvdenoise with lambda = source_map argument
%                   string - filename of the source
%   * 'fdir_out': directory to write the output (default: pwd)
%   * 'algorithm': optimization algorithm, 'quasi-newton', 'lbfgs' (default)
%   * 'num_workers': number of workers (default: 12)
%   * 'record_to': record the progress to the specified .mat file, leave '' if not recording (default: '')
% This function assumes that there is no beam loss / gain during the interaction
% Test:
% invert_shadowgraphy('H:\Projects\scripts\proton_radiography\voronoi_method\standalone\test_pic\target_map_1.png', 'num_sites', '1000', 'fdir_out', 'H:\Projects\scripts\proton_radiography\voronoi_method\standalone\test_pic', 'num_workers', '4');

function invert_shadowgraphy(filename, varargin)
    add_libs;

    totalTimeTic = tic;
    if (nargin == 0)
        error('Please specify a filename');
    elseif (mod(nargin,2) == 0)
        error('Each option must have a value accompany it.');
    end

    % set the default options
    options = {};
    options.num_sites = -1;
    options.source_map = 0;
    options.fdir_out = pwd;
    options.algorithm = 'lbfgs';
    options.num_workers = 12;
    options.record_to = '';

    % read the options specified from the user
    for (i = [1:2:length(varargin)])
        s = varargin{i}; % check the option name
        c = varargin{i+1}; % check the option value

        if (strcmp(s, 'num_sites'))
            if (ischar(c)) c = str2num(c); end;
            if (isnumeric(c))
                if (length(c) == 1)
                    options.num_sites = c;
                else
                    error('The num_sites options must be a single number');
                end
            else
                error('The num_sites options must be a single number');
            end

        elseif (strcmp(s, 'source_map'))
            if ((c == 0) | (c == '0'))
                options.source_map = 0;
            elseif (ischar(c))
                cc = str2num(c);
                if (length(cc) == 1)
                    options.source_map = cc;
                else
                    options.source_map = c;
                end
            elseif ((isnumeric(c) & (length(c) == 1)))
                options.source_map = c;
            else
                error('The source_map options must be either a string or a number');
            end

        elseif (strcmp(s, 'fdir_out'))
            if (ischar(c))
                options.fdir_out = c;
            else
                error('The option fdir_out must be a string');
            end

        elseif (strcmp(s, 'record_to'))
            if (ischar(c) & (length(c) > 3) & strcmp(c(end-3:end), '.mat'))
                options.record_to = c;
            elseif (~strcmp(c, ''))
                error('The option record_to must be a string and end in .mat');
            end

        elseif (strcmp(s, 'algorithm'))
            availAlgorithms = {'quasi-newton', 'lbfgs'};
            if (ischar(c))
                if (~string_in(c, availAlgorithms))
                    error(sprintf('The algorithm %s is not available', c));
                else
                    options.algorithm = c;
                end
            else
                error('The option algorithm must be a string');
            end

        elseif (strcmp(s, 'num_workers'))
            if (ischar(c)) c = str2num(c); end;
            if (isnumeric(c))
                if (length(c) == 1)
                    options.num_workers = c;
                else
                    error('The num_workers options must be a single number');
                end
            else
                error('The num_workers options must be a single number');
            end

        else
            error(sprintf('The options %s is not available. They are case-sensitive.', s));
        end
    end

    % parse the filename
    slashIdx = strfind(filename, '/');
    if (length(slashIdx) == 0) slashIdx = strfind(filename, '\'); end
    dotIdx = strfind(filename, '.');
    if (length(slashIdx) ~= 0) slashIdx = slashIdx(end); end;
    if (length(dotIdx) ~= 0) dotIdx = dotIdx(end); end;

    if (slashIdx == [])
        if (dotIdx == []) fname = filename;
        else fname = filename(1:dotIdx-1); end
    elseif (dotIdx == [])
        if (slashIdx == []) fname = filename;
        else fname = filename(slashIdx+1:end); end;
    else
        fname = filename(slashIdx+1:dotIdx-1);
    end

    % start the parallel pool
    disp('Creating the parallel pool workers');
    pp = gcp('nocreate'); % If no pool, do not create new one.
    if (~isempty(pp)) delete(pp); end; % if there is a pool, delete it
    cluster = parcluster('local');
    cluster.NumWorkers = options.num_workers;
    pp = parpool(cluster, options.num_workers);

    % read the file
    disp('Reading the target map file');
    img = double(imread(filename));
    img = check_image_and_convert(img);
    Npix = size(img,1) * size(img,2);

    % assign the variables according to the option values
    % num_sites
    if (options.num_sites < 0)
        options.num_sites = min(floor(0.8*Npix), 100000);
    end
    % source_map
    if (ischar(options.source_map)) % read the file
        disp('Reading the source map file');
        options.source_map = double(imread(options.source_map));
        options.source_map = check_image_and_convert(options.source_map);
    elseif (options.source_map == 0) % uniform
        disp('Assuming a uniform source map');
        options.source_map = ones(size(img));
    elseif (options.source_map > 0) % denoise from the target image
        disp(sprintf('Applying TV denoise algorithm with lambda parameter %f', options.source_map));
        options.source_map = tvdenoise(img, options.source_map, 1000);
    end
    % fdir_out
    if ((~strcmp(options.fdir_out(end), '/')) || (~strcmp(options.fdir_out(end), '\')))
        if (sum(options.fdir_out == '\') > sum(options.fdir_out == '/')) char = '\';
        else char = '/'; end;
        options.fdir_out = strcat(options.fdir_out, char);
    end

    % now let's the fun begin
    disp('Normalising the source and target maps');
    sourceMap = options.source_map / sum(options.source_map(:)) * Npix;
    targetMap = img / sum(img(:)) * Npix;

    disp('Retrieving the deflection potential from the source and target maps');
    bbbb = tic;
    [PhiI, sitesI, w] = main_inverse_extended(sourceMap, targetMap, options.num_sites, options.algorithm, 1, options.record_to);
    main_inverse_time = toc(bbbb);

    disp('Getting the reconstructed image from the retrieved potential');
    targetMapI = main_forward(sourceMap, PhiI);
    targetMapI = targetMapI/sum(targetMap(:)) * Npix;

    % save variables into a file
    disp('Saving the results into a file');
    fnameOut = strcat(options.fdir_out, fname, '-data.mat');
    save(fnameOut, 'sourceMap', 'targetMap', 'PhiI', 'targetMapI', 'main_inverse_time', 'sitesI', 'w');

    totalTime = toc(totalTimeTic);
    disp(sprintf('Finish in %d s with %d workers', totalTime, pp.NumWorkers));

    % remove the parallel pool
    delete(pp);
end

function r = string_in(s, strCells)
    r = 0;
    for (i = [1:length(strCells)])
        if (strcmp(strCells{i}, s)) r = 1; end
    end
end

function img = check_image_and_convert(img)
    if (length(size(img)) == 3) % if not grayscale, then display warning and convert it
        warning('The image is not in grayscale (1 channel). We will convert it to grayscale.');
        img = mean(img, 3);
    elseif (length(size(img)) == 4)
        error('This program cannot accept image with more than 3 dimensions');
    end
end
