%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.logL_EMCCD_lookup, high gain EMCCD noise log likelihood lookup 
% table
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

classdef logL_EMCCD_lookup < matlab.mixin.Copyable
    % Fast and reusable computation of EMCCD pixel log likelihood based on
    % lookup tables (expensive to set up, fast to evaluate).
    % This fast version uses linear interpolation in c and E, and does not
    % do derivatives.
    % 
    % constructor: obj=logL_EMCCD_lookup(EMgain,sigmaRead,cMax,Emax)
    properties (Access=public)
        % griddedInterpolant objects, and information about valid ranges
        EMgain=0;
        sigmaRead=0;         
        GIlinE=[];
    end
    methods
        function obj=logL_EMCCD_lookup(EMgain,sigmaRead,cMax,Emax)
            obj.EMgain=EMgain;
            obj.sigmaRead=sigmaRead;

            c_grid=sigmaRead*[-5:0.5:-3 -2.8:0.2:2.8 logspace(log10(3),log10(cMax/sigmaRead),200)];
            E_grid= [0 logspace(-7,log10(Emax),299)];   

            logL_EMCCD_setup(obj,EMgain,sigmaRead,c_grid,E_grid);
        end
    end
    methods (Access=protected)
        function logL_EMCCD_setup(obj,EMgain,sigmaRead,c_grid,E_grid)
            [CC,EE]=meshgrid(c_grid,E_grid);
            lnL= EMCCDfit.log_likelihood_EMCCD_brute_serial(CC,EE,EMgain,sigmaRead);            
            obj.GIlinE = griddedInterpolant({c_grid,E_grid},lnL','spline','linear');
        end
    end
    methods (Access=public)
        function y=lnL(obj,C,E)
            y=obj.GIlinE(C,E);
        end
        %function [y,dydE]=lnL_dE(obj,C,E)
        %    y=obj.GIlinE(C,E);
        %    dydE=obj.dGIlinE(C,E);
        %end
        function plotGrid(obj,figNum)
            if(exist('figNum','var'))   
                figure(figNum)
            else
                figure
            end
            clf
            hold on
            
            [x,y]=meshgrid(obj.GIlinE.GridVectors{1},obj.GIlinE.GridVectors{2}); 
            surf(x,y,obj.GIlinE.Values','edgecol','none')            
            plot3(x,y,obj.GIlinE.Values','r.','markersize',1)
            
            c=obj.GIlinE.GridVectors{1}';
            %plot3(c,0.1*ones(size(c)),obj.lnL(c,0.1*ones(size(c))),'r-','linew',1)
            c=c(c>0);
            E=c/obj.EMgain;
            plot3(c,E,obj.lnL(c,E),'r--','linew',1)
            %plot3(c,1e-1*E,obj.lnL(c,1e-1*E),'r--','linew',1)
            xlabel('c')
            ylabel('E')            
        end
    end
end
