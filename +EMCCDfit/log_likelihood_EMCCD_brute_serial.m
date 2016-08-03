function [lnL,dlnLdE]= log_likelihood_EMCCD_brute_serial(dc,E,EMgain,sigma_noise)

% This code computes the loglikelihood function for EMCCD noise model for
% the offset-subtracted data dc and photon intensity E (assuming Poisson
% photon distribution), using brute force convolution of the readout noise.
%
% ML 2015-11-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.log_likelihood_EMCCD_brute_serial, numerically compute log 
% likelihood for high gain EMCCD noise model
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


alpha = 1/EMgain;

%%% log-likelihood function
lnL = zeros(size(dc));
dlnLdE= zeros(size(dc));

dc=round(dc); % just to be sure...

% normalization constant for readout noise distribution
dcNorm=-ceil(15*sigma_noise):ceil(15*sigma_noise);
lnZdc=log(sum(exp(-0.5*(dcNorm/sigma_noise).^2)));

for ii = 1:numel(dc)

    % electron counts to consider in the numerical convolution
    if(E(ii)>0)
        Sii=max(0,floor(dc(ii)-sigma_noise*15)):max(ceil(sigma_noise*15),ceil(dc(ii)+sigma_noise*15));
        Eii=E(ii)*ones(size(Sii));
    elseif(E(ii)==0) % then no electrons are produced, and only Sii=0 is relevant
        Sii=0;
        Eii=0;
    else
       error('log likelihood not defined for E<0') 
    end
    
    % corresponding range of readout noise
    dSii=Sii-dc(ii);
    
    % log likelihoods
    lnQii=fun_lnQ(Sii,Eii,alpha);
    lnPdc=-0.5*(dSii/sigma_noise).^2;
    lnPdc=lnPdc-lnZdc;%log(sum(exp(lnPdc)));    
    % numerical convolution with noise model
    lnQP=lnQii+lnPdc;
    lnQPmax=max(lnQP);    
    lnL(ii)=lnQPmax+log(sum(exp(lnQP-lnQPmax)));

    if(nargout>1)
        % derivative of noise-free log likelihood
        dlnQii=fun_dlnQdE(Sii,Eii,alpha);
        
        % derivative of numerical convolution with noise model
        dlnQP=lnQP-lnL(ii);
        dlnQPmax=max(dlnQP);
        dlnLdE(ii)=exp(dlnQPmax)*sum(exp(dlnQP-dlnQPmax).*dlnQii);
    end
end


end
function lnQ=fun_lnQ(SS,EE,alpha)
lnQ=-inf(size(SS));
% s=0 , E>=0
lnQ(SS==0)=-EE(SS==0);
% s>0, E>0
indP=find((SS>0).*(EE>0));
xx=2*sqrt(alpha*EE(indP).*SS(indP));
% scaling: besseli(1,x)=besseli(1,x,1)*exp(abs(real(x)));
lnQ(indP) =  0.5*(log(alpha)+log(EE(indP))-log(SS(indP)));
lnQ(indP) =  lnQ(indP) -alpha*SS(indP)-EE(indP);
lnQ(indP) =  lnQ(indP) +log(besseli(1,xx,1))+abs(xx);
% the special case s>0, E=0 retains -inf

end

function dlnQ=fun_dlnQdE(SS,EE,alpha)
if(~isempty(find(EE>0,1)))
    error('log_likelihood_EMCCD cannot compute derivatives at E=0')
end
dlnQ=zeros(size(SS));
dlnQ(SS==0)=-1;

indP=find(SS>0);
xx=2*sqrt(alpha*EE(indP).*SS(indP));
dlnQ(indP) = besseli(0,xx,1)./besseli(1,xx,1).*sqrt(alpha*SS(indP)./EE(indP))-1;

end
% from http://functions.wolfram.com/Bessel-TypeFunctions/BesselI/20/ShowAll.html
% y=besseli(1,x);
% yp=besseli(0,x)-besseli(1,x)./x;
% works numerically 



























