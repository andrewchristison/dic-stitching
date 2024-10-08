function stitchImages(filepath, inputFmt, ext, options)
% STITCHIMAGES: Stitch a series of registered DIC raw images.
% 
% stitchImages(filepath, inputFmt, ext) will create a stitched output array
% from the registered images in folder filepath. Custom input formats may
% be used by changing the anonymous function inputFmt.
%
% stitchImages(___, Name, Value) is the syntax to set options.
%
% Examples:
% To stitch images in the format "C:\users\dic\e2_r3c5.png":
%   stitchImages("C:\users\dic", @(e, r, c) "e" + string(e) + "_r" +
%   string(r) + "c"+string(c), ".png.")
% To stitch images in the format "C:\users\desktop\r3c5\e0.tif"
%   stitchImages("C:\users\desktop", @(e, r, c) "r" + string(r) + "c" +
%   string(c) + "\e" + string(e), ".tif");

    arguments
        filepath string = cd + "\stitch";
        inputFmt = @(e, r, c) "e"+string(e)+"_r"+string(r)+"c"+string(c);
        ext string = ".tif";
        options.averageOverlappingRegions logical = false;
        options.outputFmt =  @(e) "e"+string(e);
        options.scaleFactor double = [1, 10];
    end
    
    try
        reg = load(fullfile(filepath, "reg.mat"));
    catch
        warning("no registration file found. calling registerImages.")
        reg = registerImages(filepath, inputFmt, ext, saveOutput = false);
    end

    idx = reg.idx; yReg = reg.y; xReg = reg.x;
    nR = size(idx.r, 1); nC =size(idx.c, 1); nE = size(idx.e, 1);

    %% MAIN
    for iE = 1:nE
        % initialize the output and its spatial reference
        stitch = nan(reg.height(iE), reg.width(iE));
        imref = imref2d([reg.height(iE), reg.width(iE)]);
        e = idx.e(iE);

        if options.averageOverlappingRegions
            oCount = nan(reg.height(iE), reg.width(iE));
        end

        % stitch all the rows and columns, starting from the bottom right
        for iRiC = (nR*nC):-1:1
            [iR, iC] = ind2sub([nR, nC], iRiC);
            r = idx.r(iR); c = idx.c(iC);

            name = fullfile(filepath, inputFmt(e, r, c) + ext);
            image = rescale(imread(name), 0, 1);
            [rS] = [1, size(image, 1)] + yReg(iR, iC, iE);
            [cS] = [1, size(image, 2)] + xReg(iR, iC, iE);
            [rS, cS] = worldToSubscript(imref, cS, rS);

            if options.averageOverlappingRegions
                oCount(rS(1):rS(2), cS(1):cS(2)) = ...
                    sum(cat(3, oCount(rS(1):rS(2), cS(1):cS(2)), ...
                    ones(size(image))), 3, "omitnan");
                stitch(rS(1):rS(2), cS(1):cS(2)) = ...
                    sum(cat(3, stitch(rS(1):rS(2), cS(1):cS(2)), ...
                    image), 3, "omitnan");
            else
                stitch(rS(1):rS(2), cS(1):cS(2)) = image;
            end
        end

        if options.averageOverlappingRegions
            stitch = stitch ./ oCount;
        end

        % save the output images for this strain level
        sF = options.scaleFactor(:);
        for s = 1:size(sF, 1)
            if sF(s) == 1
                fname = options.outputFmt(iE);
                file = fullfile(filepath, fname + ext);
                imwrite(stitch, file);
            else
                fname = options.outputFmt(iE + "_rescale" + string(s));
                scale = imresize(stitch, 1/sF(s), method = "lanczos2");
                file = fullfile(filepath, fname + ext);
                imwrite(scale, file);
            end
        end
    end
end