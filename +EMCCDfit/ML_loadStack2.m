function [rawData, stackSize] = ML_loadStack2(stackName,maxFrames)
% Reads in a tif-stack from stackName. If specified, only first maxFrames
% frames are read. This function is based on code by Fredrik Persson.
%
% ML 2014-06-19

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.ML_loadStack2.m, load tif-stacks
% =========================================================================
% 
% Copyright (C) 2016 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by linking or combining it
%  with Matlab or any Matlab toolbox, the licensors of this Program grant you 
%  additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code


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
