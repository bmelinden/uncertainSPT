% switch x and y axes in the active plot (gca). Lines, errorbars, text,
% axes, and labels are affected.
%
% ML 2016-10-27

% data
ch=get(gca,'children');
for m=1:numel(ch)
    switch class(ch(m))
        case 'matlab.graphics.chart.primitive.Line'
            x=get(ch(m),'xdata');
            y=get(ch(m),'ydata');
            set(ch(m),'xdata',y,'ydata',x);
            disp([class(ch(m)) ' :  switched.'] )
        case 'matlab.graphics.chart.primitive.ErrorBar'
            x=get(ch(m),'xdata');
            y=get(ch(m),'ydata');
            dxm=get(ch(m),'XNegativeDelta');
            dxp=get(ch(m),'XPositiveDelta');
            dym=get(ch(m),'YNegativeDelta');
            dyp=get(ch(m),'YPositiveDelta');
            
            set(ch(m),'xdata',y,'XNegativeDelta',dym,'XPositiveDelta',dyp,...
                      'ydata',x,'YNegativeDelta',dxm,'YPositiveDelta',dxp);
   
            disp([class(ch(m)) ' :  switched.'] )
        case 'matlab.graphics.primitive.Text'
            xyzp=get(ch(m),'position');
            v  =get(ch(m),'rotation');
            
            % switch
            xyzp=xyzp([2 1 3]);
            v=0*(v==90)+90*(v==0);            
            set(ch(m),'position',xyzp,'rotation',v);
            
            disp([class(ch(m)) ' :  switched.'] )
        otherwise
            disp([class(ch(m)) ' :  unchanged.'] )                
    end
end
    
xprop={'xlabel','xlim','xscale'};%,'xtick','xticklabel'};
yprop={'ylabel','ylim','yscale'};%,'ytick','yticklabel'};
for m=1:numel(xprop)
   x=get(gca,xprop{m});
   y=get(gca,yprop{m});
   set(gca,xprop{m},y,yprop{m},x)   
end



