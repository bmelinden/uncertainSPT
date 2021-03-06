function displayParameter(Pmle,ncd,Dscale,ro)
% EMhmm.displayParameter(Pmle,ncd,Dscale,ro)
%
% Displays EMhmm model parameters.
% 
% Pmle  : parameter struct to display as central estimate, see
%         EMhmm.parameterEstimate. 
% ncd   : number of characters and decimals in numeric strings. 
%         ncd(1) > ncd(2) is strongly recommended. Default : ncd = [6 3],
%         i.e., 6 characters, 3 decimal places.
% Dscale: rescale units of step variance (e.g., Dscale=1e-6 to convert 
%         nm^2/s to um^2/s), affects lambda and all variables ending in D.
%         Default 1.
% ro    : optional true/false indicator for restricted output of only D,
%         pOcc, and mean dwell times (default: false)
%
% ML 2017-01-27
 
%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displayParameter, display parameters of diffusive HMM
% =========================================================================
% 
% Copyright (C) 2017 Martin Lindén
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

%% start actual code
if( ~exist('ncd','var') || isempty(ncd) )
    ncd=[6 3];
end
if( ~exist('Dscale','var') || isempty(Dscale) )
    Dscale=1;
end

f=fieldnames(Pmle);
% check for restricted output
if( exist('ro','var') && ~isempty(ro) && ro)
    f=f(ismember(f,{'D','pOcc','dwellTime'}));
end

varLength=1;
for k=1:length(f)
    varLength=max(varLength,length(f{k}));
end
disp('parameter : ')
floatString=[' %' int2str(ncd(1)) '.' int2str(ncd(2)) 'f | '];
for k=1:length(f)
    P=Pmle.(f{k});
    for r=1:size(P,1)
        fprintf([' %' int2str(varLength) 's : '],f{k})
        for c=1:size(P,2)
            if(strcmp(f{k}(end),'D') || strcmp(f{k},'lambda'))
                fprintf(floatString,P(r,c)*Dscale)
            else
                fprintf(floatString,P(r,c))
            end
        end
        fprintf('\n')
    end
end
