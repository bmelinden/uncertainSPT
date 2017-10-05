function [rawData, stackSize] = ML_loadStack2(stackName,maxFrames)
% Reads in a tif-stack from stackName. If specified, only first maxFrames
% frames are read. This function is based on code by Fredrik Persson.
%
% ML 2014-06-19

[sPath,sName,sExt]=fileparts(stackName);
if( strfind(stackName,'.tif') )
    
    % Calc. new variables
    imInfo = imfinfo(stackName);
    imWidth = imInfo(1).Width;
    imHeight = imInfo(1).Height;
    stackSize = numel(imInfo);
    
    if(exist('maxFrames','var') && ~isempty(maxFrames) && maxFrames>0)
        stackSize=min(stackSize,maxFrames);
    end
    
    rawData = zeros(imHeight, imWidth, stackSize);
    
    t = Tiff(stackName, 'r');
    for inFrame = 1:stackSize
        inFrame;
        rawData(:, :, inFrame) = t.read();
        if inFrame<stackSize
            t.nextDirectory();
        end
    end
    t.close();
    clear 'im*' 't' 'inFrame' 'filename' 'pathname';  
else
    error('ML_loadStack2.m can only handle tif stacks yet.')
end
