clear;
proportion='0';
%%
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
Force=dlmread([FRACTION 'force_info.asc'],' ');Force=Force([1:10:2400,2400],:);
Disp=dlmread([FRACTION 'displacement_info.asc'],' ');Disp=Disp(1:241,:);
pos_vec=load('../../Simulation_DATA/DirectShearTest/ball_pos_PFC/Ball_pos.txt','r');
info=textread('../../Simulation_DATA/DirectShearTest/ball_pos_PFC/Ball&Wall_info.txt','%f');
x_min=info(1);
x_max=info(2);
y_min=info(3);
y_max=info(4);
rhi=info(6);
x_min=-6*ceil(-x_min/(6*rhi))*rhi;
x_max=6*ceil(x_max/(6*rhi))*rhi;
%%
ShearStress=0.5*(Force(:,3)-Force(:,5))/(x_max-x_min);
NormalStress=0.5*(Force(:,6)-Force(:,4))/(x_max-x_min);
StressRatio=ShearStress./NormalStress;
VerticalDisp=(Disp(:,6)-Disp(1,6))-(Disp(:,4)-Disp(1,4));
VerticalStrain=VerticalDisp./(y_max-y_min);
SolidVol=sum(pi*pos_vec(:,3).^2);
BoxVol=(x_max-x_min)*(y_max-y_min+Disp(:,6)+VerticalDisp+0.625*rhi);
SolidFraction=(BoxVol-SolidVol)./SolidVol;
[ShearStrain]=Shearstrain(proportion);

%%
plot(ShearStrain,VerticalStrain);hold on;