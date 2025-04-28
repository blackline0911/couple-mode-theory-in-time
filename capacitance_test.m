Na = 3.31e+17;
Nd =  6.44e+17 ;

eps_Si = 3.47^2;
eps_clad = 1.44^2;
eps_box = eps_clad;
e0 = 8.854e-14;
q = 1.6e-19;
ni = 1.07e10;
La = 2*pi*5*1e-4;
hrib = 0.2114*1e-4;


V0 = 0.0259*log(Na*Nd/ni^2);
V = linspace(-2,0,1000);
W = sqrt( 2*e0*eps_Si*(V0-V) / (q*Na*Nd*(Nd+Na)) ) * (Na+Nd);
Wp = W*Nd/(Na+Nd);
Wn = W*Na/(Na+Nd);
tpb = 0.25e-4;
tnb = 0.25e-4;
C_parallel = e0*eps_Si*hrib*La./W;
Cf = La*e0*(eps_clad+eps_box)/2/pi*log(2*pi*hrib./W);
k = sqrt( W.*(W+tpb).*(W+tnb+tpb)./ (W+tnb) );
k_pround = sqrt(1-k.^2);
C_top = La*e0*eps_clad*ellipke(k_pround) ./ellipke(k);
C_bottom = C_top;
Cj = C_parallel+Cf+C_top+C_bottom; 
Cj_pdk = [23.6*1e-15,20*1e-15,18.5*1e-15];

figure('Color','w')
plot([0,-1,-2],Cj_pdk,'b.','Marker','o')
hold on
plot(V,Cj,"Linewidth",2);
hold on
% plot(V,C_parallel,"Linewidth",2)
xlabel("Voltage");
title("Junction capacitance");
legend("PDK data","Cj deduced");
