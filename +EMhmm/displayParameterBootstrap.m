function displayParameterBootstrap(Pmle,BS,ncd)
% EMhmm.displayParameterBootstrap(Pmle,BS,ncd)
%
% Displays parameter +- bootstrap std. err.
% 
% Pmle  : parameter struct, see EMhmm.parameterEstimate
% BS    : bootstrap parameter struct, see EMhmm.parameterBootstrap
% ncd   : number of characters and decimals in numeric strings. 
%         ncd(1) > ncd(2) is recommended. efault : ncd = [6 3], i.e., 
%         6 characters, 3 decimal places.
%
% ML 2016-08-19
 
%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displayParameterBootstrap, display parameters of diffusive HMM
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


if( ~exist('ncd','var') || isempty(ncd) )
    ncd=[6 3];
end

f=fieldnames(Pmle);
varLength=1;
for k=1:length(f)
    varLength=max(varLength,length(f{k}));
end
disp('MLE parameters +- std err : ')
floatString=[' %' int2str(ncd(1)) '.' int2str(ncd(2)) 'f +- %' int2str(ncd(1)) '.' int2str(ncd(2)) 'f | '];
for k=1:length(f)
    P=Pmle.(f{k});
    dP=std(BS.(f{k}),[],3);
    for r=1:size(P,1)
        fprintf([' %' int2str(varLength) 's : '],f{k})
        for c=1:size(P,2)
            fprintf(floatString,P(r,c),dP(r,c))
        end
        fprintf('\n')
    end
end
