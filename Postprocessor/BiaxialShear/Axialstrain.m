function [weyy]=Axialstrain(proportion)
%%
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
Disp=dlmread([FRACTION 'displacement_info.asc'],' ');
info=textread('../../Simulation_DATA/BiaxialShearTest/ball_pos_PFC/Ball&Wall_info.txt','%f');
y_min=info(3);
y_max=info(4);

%%
ly_0=y_max-y_min-Disp(1,4)+Disp(1,6);
wly=ly_0-(Disp(:,4)-Disp(1,4))+(Disp(:,6)-Disp(1,6));
weyy=-(wly-ly_0)./ly_0;
end