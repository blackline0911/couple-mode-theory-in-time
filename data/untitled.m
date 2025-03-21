close all
%% parameter setting
%% Note that in coupled mode theory in time(CMT) , we must ensure decay rate is far smaller than resonant frequency w0.
%% And the CMT is valid when decay rate is small enough that the spectral broadening does not overlap the near spectral mode.  
gamma = 0.9581231653602788;
kappa_1 = sqrt(1-gamma.^2);
a = 0.9867369707049318;
r = 50e-6;
ng =3.93;
neff = [2.51464,2.5147];
c = 299792458;
tu_e1 = -2*pi*r.*ng./(c*log(sqrt(1-kappa_1^2)));
tu_o = -2*pi*r.*ng./(c*log(a));  % photon life time
wl = linspace(1.32079,1.32086,10000)*1e-6;
w = 2*pi*c./wl;
% w0 = 102*pi*c./(neff*2*pi*r)
decay_rate = 1/tu_e1+1/tu_o
lambda0 = 1.32082;
w0 = 2*pi*c./lambda0/1e-6;

data = readmatrix("T.txt");
lambda = readmatrix("wl.txt");
%% coupling equation in time

% T_time = ( 1/tu_e1-1/tu_o-1i*(w-w0) )./( 1/tu_e1+1/tu_o+1i*(w-w0) );
figure('color','w')
for i=1:1
    T_time = ( -1/tu_e1+1/tu_o+1i*(w-w0(i)) )./( 1/tu_e1+1/tu_o+1i*(w-w0(i)) );
    % t = sqrt(1-kappa_1^2);
    % T_time = ( 3/tu_e1+1/tu_o+1i*(w-w0) )./(1/tu_e1+1/tu_o+1i*(w-w0));
    
    plot(wl,abs(T_time).^2,'Linewidth',2);
    [min_value, min_idx] = min(abs(T_time).^2);
    x_min_1 = wl(min_idx) % 對應的 X 軸值
    hold on
end
[min_value_data, min_idx] = min(data);
x_min = lambda(min_idx) % 對應的 X 軸值
plot(lambda*1e-6,data,'Linewidth',2)
x_min_1*1e6- x_min
min_value-min_value_data
hold off

% 
data_1=readmatrix("T.txt");
data_2=readmatrix("T_NL.txt");
plot(lambda,10*log10(data_1),'Linewidth',2)
hold on
plot(lambda,10*log10(data_2),'Linewidth',2)