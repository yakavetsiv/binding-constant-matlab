function [k,go,out,xx,curve]=FitModel11(x,y)
%задание граничный условий и параметров нелинейной регрессии
fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On',...
    'Lower',[1],'Upper',[Inf],'DiffMinChange',1e-03,...
    'TolFun',1e-012,'TolX',1e-11,'MaxIter',1e12);
ok_ = isfinite(x) & isfinite(y);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs',...
        'Ignoring NaNs and Infs in data.' );
end
%начальные условия
st_ = [6e6];
p0=0;
%задание аппроксимационной функции
set(fo_,'Startpoint',st_);
ft_ = fittype(@(k1,p0,x)ff11(k1,p0,x),...
    'problem',{'p0'},'dependent',{'y'},'independent',{'x'});

%аппроксимация
[cf_,gof] = fit(x(ok_),y(ok_),ft_,fo_,'problem',{p0});
%полученные параметры
par_=coeffvalues(cf_);
k=[par_(1)];
go=[gof.adjrsquare,gof.rmse];
ci = confint(cf_,0.95);
out=[par_(1);gof.adjrsquare;gof.rmse;ci];

%построение графика экспериментальных данных и полученной кривой
h=semilogx(x,y);
set(h,'Marker','.');
axis([1e-9 1e-3 0 1]);
hold on;
xx=logspace(-8,-4);
curve=ff11(k(1),p0,xx);
plot(xx,curve,'-r');
hold off;

%модельная функция комплексообразования
function y=ff11(k1,p0,xdata)
p1=1;
G0=xdata;
y1=p0+p1*k1*G0;
y2=1+k1*G0;
y=y1./y2;
end
end
