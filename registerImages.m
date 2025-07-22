function reg = registerImages(filepath, inputFmt, ext, options)
% REGISTER IMAGES: registers a series of DIC images with cross correlation
%
% reg = registerImages(filepath, inputFmt, ext) will automatically find all
% files with the specified format in the filepath.
%
% Examples:
% To save registrations in a .mat file:
%   registerImages(____, saveOutput = true);
% To register images in the format "C:\users\dic\e2_r3c5.png":
%   registerImages("C:\users\dic", @(e, r, c) "e" + string(e) + "_r" +
%   string(r) + "c"+string(c), ".png.")

    arguments
        filepath string = cd + "\stitch";
        inputFmt = @(e, r, c) "e"+string(e)+"_r"+string(r)+"c"+string(c);
        ext string = ".tif";
        options.saveOutput logical = true;
        options.displayProgress logical = true;
        options.cropTemplate double = 1024;
    end

    [names, idx] = parseFilenames(filepath, inputFmt, ext);
    reg.idx = idx;              [nR, nC, nE] = size(names);
    xGlobal = nan(nR, nC, nE);  yGlobal = nan(nR, nC, nE);
    h = zeros(1, nE);           w = zeros(1, nE);

    cropTemplate = options.cropTemplate;

    %% MAIN
    for iE = 1:nE

        % load and register all the columns
        if options.displayProgress
            fprintf("Strain level %d. Cols registered:", idx.e(iE));
        end
        yC = nan(nR, nC); xC = nan(nR, nC); hC = zeros(1, nC);
        for iC = 1:nC
            [yC(:, iC), xC(:, iC), hC(iC), ~] = ...
                register(names(:, iC, iE), cropTemplate);
            if options.displayProgress
                fprintf(" %d", idx.c(iC));
            end
        end
        
        % load and register all the rows
        if options.displayProgress
            fprintf("\nStrain level %d. Rows registered:", idx.e(iE));
        end
        yR = nan(nR, nC); xR = nan(nR, nC); wR = zeros(nR, 1);
        for iR = 1:nR
            [yR(iR, :), xR(iR, :), ~, wR(iR)] = ...
                register(names(iR, :, iE), cropTemplate);
            if options.displayProgress
                fprintf(" %d", idx.r(iR));
            end
        end
        if options.displayProgress
            fprintf("\n");
        end

        yGlobal(:, :, iE) = mean(cat(3, yR + yC(:, 1), yC + yR(1, :)), 3);
        xGlobal(:, :, iE) = mean(cat(3, xR + xC(:, 1), xC + xR(1, :)), 3);
        h(iE) = ceil(max(yGlobal(end, :, iE) + hC) + 1);
        w(iE) = ceil(max(xGlobal(:, end, iE) + wR) + 1);
    end

    if options.displayProgress
        fprintf("Registration complete. Saving to file.\n");
    end

    reg.y = yGlobal; reg.x = xGlobal; reg.height = h; reg.width = w;

    if options.saveOutput == true
        outputFile = fullfile(filepath, "reg.mat");
        save(outputFile, "-struct", "reg");
    end
end

% implements registration for a row or column of images
function [yPos, xPos, h, w] = register(names, cropTemplate)

    % initialize loop variables
    i = 1;
    template = imread(names(1));
    numFiles = numel(names);

    % initialize output arrays
    xLocal = nan(numFiles, 1); xLocal(1) = 0;
    yLocal = nan(numFiles, 1); yLocal(1) = 0;
    xPos = xLocal; yPos = yLocal;

    while(i < numFiles)
        image = template;
        template = imread(names(i+1));
        
        % local registration for a pair of images
        ct = imcrop(template, [1, 1, cropTemplate, cropTemplate]);
        cc = normxcorr2(ct, image);
        [yPeak,xPeak] = find(cc==max(cc(:)));
        yPeak = yPeak - size(ct, 1);
        xPeak = xPeak - size(ct, 2);

        % global registration
        yLocal(i+1) = yPeak; yPos(i+1) = yPeak + yPos(i);
        xLocal(i+1) = xPeak; xPos(i+1) = xPeak + xPos(i);

        i = i + 1;
    end

    % move coordinates to have a (0, 0) origin
    y0 = min(yPos); yPos = yPos - y0;
    x0 = min(xPos); xPos = xPos - x0;
    h = size(image, 1); w = size(image, 2);
end

function [names, idx] = parseFilenames(filepath, inputFmt, ext)
    % parse all the possible files
    ls = dir(filepath + "\**\" + "*"+ext);
    ls = struct2table(ls);
    testNames = fullfile(string(ls.folder), string(ls.name));

    % parse and tokenize all files to be stitched in the folder
    expression = regexptranslate("escape", filepath + "\") ...
               + regexptranslate("escape", inputFmt("<e>","<r>","<c>"))...
               + regexptranslate("wildcard", ext);
    expression = replace(expression, ["<e>", "<r>", "<c>"], ...
                   ["(?<e>\d*)", "(?<r>\d*)", "(?<c>\d*)"]);
    [names, tokens] = regexp(testNames, expression, 'match', 'names');
    names = [names{:}]; names = names(:);

    % sort all the tokens
    tokens = convertvars(struct2table([tokens{:}]), @isstring, "double");
    [tokens, sorted] = sortrows(tokens, [1, 3]);
        
    % get all the unique values of e, r, c.
    idx.e = unique(tokens.e);
    idx.r = unique(tokens.r);
    idx.c = unique(tokens.c);

    % reshape the array of file names to be 3D
    nE = size(idx.e, 1); nR = size(idx.r, 1); nC = size(idx.c, 1);
    names = reshape(names(sorted), [nR, nC, nE]);
end