clear;
proportion='30';
%%
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
Force=dlmread([FRACTION 'force_info.asc'],' ');Force=Force([1:10:350,350],:);
Disp=dlmread([FRACTION 'displacement_info.asc'],' ');Disp=Disp([1:10:350,350],:);
info=textread('../../Simulation_DATA/BiaxialShearTest/ball_pos_PFC/Ball&Wall_info.txt','%f');
x_min=info(1);
x_max=info(2);
y_min=info(3);
y_max=info(4);
%%
lx_0=x_max-x_min-Disp(1,7)+Disp(1,9);
ly_0=y_max-y_min-Disp(1,4)+Disp(1,6);
wlx=lx_0-(Disp(:,7)-Disp(1,7))+(Disp(:,9)-Disp(1,9));
wly=ly_0-(Disp(:,4)-Disp(1,4))+(Disp(:,6)-Disp(1,6));
wsxx=(Force(:,9)-Force(:,7))./(2.*wly);
wsyy=(Force(:,6)-Force(:,4))./(2.*wlx);
p=0.5*(wsxx+wsyy);
q=0.5*(wsyy-wsxx);
wexx=(wlx-lx_0)./lx_0;
weyy=-(wly-ly_0)./ly_0;
wevol=(wlx.*wly-lx_0.*ly_0)./(lx_0.*ly_0);

%%
plot(weyy,q./p);hold on;