function stitchDic(filepath, inputFmt, ext, options)
% STITCHDIC: Stitch a series of registered, correlated DIC fields.
%
% STITCHDIC(filepath, inputFmt, ext, options) will create a stitched output
% array from the registered .mat images in folder filepath. Custom input
% formats may be used by changing the anonymous function inputFmt.
%
% Data is assumed to be a Matlab struct with displacement and strain data
% as fields e.g., dic.exx, dic.u, dic.x
%
% stitchDic(___, Name, Value) is the syntax to set options.
% 
% All DIC fields must contain a correlation flag as set in the options. By
% default, this flag is set to the sigma field from VIC-2D output.
%
% Examples:
% To stitch DIC data in the format "C:\users\desktop\r3c5\e1.mat":
%   stitchDic("C:\users\desktop", @(e, r, c) "r" + string(r) + "c" +
%   string(c) + "\e" + string(e), ".mat");
% To change the displacement field names to "U" and "V":
%   stitchDic(___, displacementFields = {"U, V"});
% To output data in the deformed configuration:
%   stitchDic(___, outputDeformedConfig = true);

    arguments
        filepath string = cd + "\stitch";
        inputFmt = @(e, r, c) "e"+string(e)+"_r"+string(r)+"c"+string(c);
        ext string = ".mat";
        options.corrFlag string = "sigma";
        options.corrFail = -1;
        options.outputFmt =  @(e) "e"+string(e)+"_stitched";
        options.displacementFields = {"u", "v"};
        options.deformationFields = {"exx", "eyy", "exy"};
        options.smoothFields = {"u", "v", "exx"};
        options.defineReferenceConfig double = 0;
        options.outputDeformedConfig logical = false;
        options.filterConfidence logical = false;
    end

    try
        reg = load(fullfile(filepath, "reg.mat"));
    catch
        warning("no registration file found. calling registerImages.")
        reg = registerImages(filepath, inputFmt, ext, saveOutput = false);
    end

    % get the reference/deformed registrations
    idx = reg.idx; e0 = (idx.e == options.defineReferenceConfig);
    idx.e = idx.e(~e0);
    xReg_0 = reg.x(:, :, e0); yReg_0 = reg.y(:, :, e0);
    reg.x = reg.x(:, :, ~e0); reg.y = reg.y(:, :, ~e0);

    % get the step size and struct fields from the first input file
    fields = [options.deformationFields(:); 
              options.displacementFields(:); 
              options.corrFlag];
    u = options.displacementFields{1}; v = options.displacementFields{2};
    nF = size(fields(:), 1);

    nR = size(idx.r, 1); nC = size(idx.c, 1); nE = size(idx.e, 1);

    %% MAIN
    for iE = 1:nE
        e = idx.e(iE);
        xReg_E = reg.x(:, :, iE); yReg_E = reg.y(:, :, iE);

        % load the (0, 0)/origin DIC file
        name = fullfile(filepath,inputFmt(e, idx.r(1), idx.c(1)) + ext);
        corrOrigin = load(name);

        % initialize the output and spatial reference
        stepSize = corrOrigin.x(1, 2) - corrOrigin.x(1, 1);
        [X, Y] = meshgrid(1:stepSize:reg.width(e0), 1:stepSize:reg.height(e0));
        imref = imref2d(size(X), stepSize, stepSize);
        for f = 1:nF
            stitch.(fields{f}) = nan(imref.ImageSize);
        end

        % corrGlobal should flag if correlation has failed
        corrGlobal = options.corrFail * ones(imref.ImageSize);
        for iRiC = (nR*nC):-1:1
            [iR, iC] = ind2sub([nR, nC], iRiC);
            r = idx.r(iR); c = idx.c(iC);
            name = fullfile(filepath, inputFmt(e, r, c) + ext);
            try
                corrLocal = load(name);

                % adjust displacements
                corrLocal.(u) = corrLocal.(u)...
                            + (xReg_E(iR, iC) - xReg_0(iR, iC));
                corrLocal.(v) = corrLocal.(v)...
                            + (yReg_E(iR, iC) - yReg_0(iR, iC));

                % shift corrLocal into global coordinates
                corrLocal.x = corrLocal.x + xReg_0(iR, iC);
                corrLocal.y = corrLocal.y + yReg_0(iR, iC);

                % filter for the best correlating subsets
                [maskL, maskG] = ...
                    filterSubsets(corrLocal, corrGlobal, imref, options);
                corrGlobal(maskG) = corrLocal.(options.corrFlag)(maskL);

                % stitch all the fields
                for f = 1:nF
                    stitch.(fields{f})(maskG) = corrLocal.(fields{f})(maskL);
                end
            catch
                warning("%s failed to stitch.", name);
            end
        end

        sF = options.smoothFields;
        for i = 1:size(sF(:), 1)
            stitch.(sF{i}) = smoothdata(stitch.(sF{i}), 1, "gaussian", 5);
            stitch.(sF{i}) = smoothdata(stitch.(sF{i}), 2, "gaussian", 5);
        end

        % NOTE: interpolating missing displacement data is very slow
        % for deformed config, smoothing displacements is recommended
        if options.outputDeformedConfig
            u = fillmissing2(stitch.(u), "linear", ...
                MissingLocations = isnan(stitch.(u)));
            v = fillmissing2(stitch.(v), "linear", ...
                MissingLocations = isnan(stitch.(v)));
            D = cat(3, u, v) ./ stepSize; D(isnan(D)) = 0;
            for f = 1:nF
                    stitch.(fields{f}) = imwarp(stitch.(fields{f}), D, ...
                                         "nearest", FillValues= nan);
            end
        end
        stitch.x = X; stitch.y = Y;

        % save outputs for this strain level
        %file = fullfile(filepath, options.outputFmt(iE) + ext);
        %save(stitch, file);
    end
end

function [maskL, maskG] = filterSubsets(corrLocal, corrGlobal, imref, options)
    flag = options.corrFlag;
    fail = options.corrFail;

    % reference corrLocal to the global coordinate system
    xG = [corrLocal.x(1, 1), corrLocal.x(end, end)];
    yG = [corrLocal.y(1, 1), corrLocal.y(end, end)];
    [rS, cS] = worldToSubscript(imref, xG, yG);

    corrGlobal_iRiC = fail * ones(size(corrGlobal));
    corrGlobal_iRiC(rS(1):rS(2), cS(1):cS(2)) = corrLocal.(flag);

    % create masks to merge iRiC into corrGlobal
    % maskG: subsets correlated in iRiC
    % maskC: subsets where iRiC had better confidence than corrGlobal
    maskG = (corrGlobal_iRiC ~= fail);
    if options.filterConfidence
        maskC = (corrGlobal == fail) | (corrGlobal_iRiC < corrGlobal);
        maskG = maskG & maskC;
    end
    maskL = maskG(rS(1):rS(2), cS(1):cS(2));
end