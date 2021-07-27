%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invert the grayscale shadowgraphy image
% Input:
%   * filename: string to a filename
% Options:
%   * 'source_map': 0 - uniform with the same size as the file (default)
%                   >0 - using tvdenoise with lambda = source_map argument
%                   string - filename of the source
%   * 'fdir_out': directory to write the output (default: pwd)
%   * 'record_to': record the progress to the specified .mat file, leave '' if
%                  not recording (default: '')
% This function assumes that there is no beam loss / gain during the interaction
% Test:
% fast_invert_shadowgraphy('test-figures/test-100pix.png', ...
%                          'fdir_out', '/your/intended/directory');
%
% Author: Muhammad F. Kasim (University of Oxford, 2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fast_invert_shadowgraphy(filename, varargin)
    add_libs;

    totalTimeTic = tic;
    if (nargin == 0)
        error('Please specify a filename');
    elseif (mod(nargin,2) == 0)
        error('Each option must have a value accompany it.');
    end

    % set the default options
    options = {};
    options.source_map = 0;
    options.fdir_out = pwd;
    options.record_to = '';

    % read the options specified from the user
    for (i = [1:2:length(varargin)])
        s = varargin{i}; % check the option name
        c = varargin{i+1}; % check the option value

        if (strcmp(s, 'source_map'))
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

    % read the file
    fprintf('Reading the target map file\n');
    img = double(imread(filename));
    img = check_image_and_convert(img);
    Npix = size(img,1) * size(img,2);

    % assign the variables according to the option values
    % source_map
    if (ischar(options.source_map)) % read the file
        fprintf('Reading the source map file\n');
        options.source_map = double(imread(options.source_map));
        options.source_map = check_image_and_convert(options.source_map);
    elseif (options.source_map == 0) % uniform
        fprintf('Assuming a uniform source map\n');
        options.source_map = ones(size(img));
    elseif (options.source_map > 0) % denoise from the target image
        fprintf('Applying TV denoise algorithm with lambda parameter %f\n', options.source_map);
        options.source_map = tvdenoise(img, options.source_map, 1000);
    end
    % fdir_out
    if ((~strcmp(options.fdir_out(end), '/')) || (~strcmp(options.fdir_out(end), '\')))
        if (sum(options.fdir_out == '\') > sum(options.fdir_out == '/')) char = '\';
        else char = '/'; end;
        options.fdir_out = strcat(options.fdir_out, char);
    end

    % now let's the fun begin
    fprintf('Normalising the source and target maps\n');
    sourceMap = options.source_map / sum(options.source_map(:)) * Npix;
    targetMap = img / sum(img(:)) * Npix;

    fprintf('Retrieving the deflection potential from the source and target maps\n');
    bbbb = tic;
    PhiI = main_inverse(sourceMap, targetMap,options);
    main_inverse_time = toc(bbbb);

    fprintf('Getting the reconstructed image from the retrieved potential\n');
    targetMapI = main_forward(sourceMap, PhiI);
    targetMapI = targetMapI/sum(targetMap(:)) * Npix;

    % save variables into a file
    fprintf('Saving the results into a file\n');
    fnameOut = strcat(options.fdir_out, fname, '-data.mat');
    save(fnameOut, 'sourceMap', 'targetMap', 'PhiI', 'targetMapI', 'main_inverse_time');

    subplot(2,2,1);
    imagesc(sourceMap);
    title('Source');
    colorbar;
    subplot(2,2,2);
    imagesc(targetMap);
    title('Original target');
    colorbar;
    subplot(2,2,3);
    imagesc(PhiI);
    title('Deflection potential');
    colorbar;
    subplot(2,2,4);
    imagesc(targetMapI);
    title('Reconstructed target');
    colorbar;

    totalTime = toc(totalTimeTic);
    fprintf('Finish in %ds\n', totalTime);
end

function r = string_in(s, strCells)
    r = 0;
    for (i = [1:length(strCells)])
        if (strcmp(strCells{i}, s)) r = 1; end
    end
end

function img = check_image_and_convert(img)
    if (length(size(img)) == 3) % if not grayscale, then fprintflay warning and convert it
        warning('The image is not in grayscale (1 channel). We will convert it to grayscale.');
        img = mean(img, 3);
    elseif (length(size(img)) == 4)
        error('This program cannot accept image with more than 3 dimensions');
    end
end
