function [ShearStrain]=Shearstrain(proportion)
%%
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
Disp=dlmread([FRACTION 'displacement_info.asc'],' ');Disp=Disp(1:241,:);
info=textread('../../Simulation_DATA/DirectShearTest/ball_pos_PFC/Ball&Wall_info.txt','%f');
y_min=info(3);
y_max=info(4);

%%
HorizontalDisp=(Disp(:,5)-Disp(1,5))-(Disp(:,3)-Disp(1,3));
VerticalDisp=(Disp(:,6)-Disp(1,6))-(Disp(:,4)-Disp(1,4));
ShearStrain=HorizontalDisp./(y_max-y_min-VerticalDisp);
end