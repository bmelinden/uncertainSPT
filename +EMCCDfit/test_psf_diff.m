% numerical test of partial derivatives for an asymmetric Gaussian psf
% model
clear

dp=linspace(-1,1,2000);
p0={randn,randn,-1.5,1,2,1,pi*randn };
x=randn;
y=1+randn;

%%
for k=1:7
    pk=p0;
    g=p0;
    for m=1:7
        pk{m}=pk{m}+zeros(size(dp));
        g{m} =zeros(size(dp));
    end
    pk{k}=pk{k}+dp;
    for m=1:length(dp)
        [E(m),g{1}(m),g{2}(m),g{3}(m),g{4}(m),g{5}(m),g{6}(m),g{7}(m)]=....
            EMCCDfit.psf_diff_asymgauss_angle(x,y,...
        pk{1}(m),pk{2}(m),pk{3}(m),pk{4}(m),pk{5}(m),pk{6}(m),pk{7}(m));
    end
    figure(1)
    clf
    hold on
    P=pk{k};
    dEdp=diff(E)./diff(P);
    
    plot(P,g{k},'k')
    plot(P(1:end-1)+0.5*mean(diff(dp)),dEdp,'r')
    pause(0.1)
    
end