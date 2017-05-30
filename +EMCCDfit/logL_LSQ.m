%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.logL_LSQ, least-square likelihood function
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

classdef logL_LSQ < matlab.mixin.Copyable
    % Fast and reusable computation of EMCCD pixel log likelihood based on
    % lookup tables (expensive to set up, fast to evaluate).
    % This fast version uses linear interpolation in c and E, and does not
    % do derivatives.
    % 
    % constructor: obj=logL_EMCCD_lookup(EMgain,sigmaRead,cMax,Emax)
    properties (Access=public)
    end
    methods
        function obj=logL_LSQ()
        end
    end
    methods (Access=public)
        function y=lnL(obj,C,E)
            y=-(C-E).^2;
        end
    end
end
