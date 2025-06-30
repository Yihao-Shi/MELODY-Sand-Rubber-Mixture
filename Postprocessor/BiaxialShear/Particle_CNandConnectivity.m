clear
proportion='10';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
contact=dlmread([FRACTION 'contacting_info.asc'],' ');
contact=contact(:,3:end-1);

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

%% Coordination Number and Connectivity
Z_soft=zeros(size(contact,1),1);
Z_rigid=zeros(size(contact,1),1);
Zsoft_m=zeros(size(contact,1),1);
Zrigid_m=zeros(size(contact,1),1);
Z=zeros(size(contact,1),1);
Z_m=zeros(size(contact,1),1);

for timestep=1:160%size(contact,1)
    %Connectivity per sec
    fai_soft=contact(timestep,1:soft_num).';
    fai_rigid=contact(timestep,soft_num+1:particle_num).';
    
    %average coordinate
    C_soft=sum(fai_soft);
    Z_soft(timestep)=C_soft/soft_num;
    C_rigid=sum(fai_rigid);
    Z_rigid(timestep)=C_rigid/rigid_num;
    Z(timestep)=mean(contact(timestep,:));
    
    %mechanical coordinate
    Nsoft_0=length(find(fai_soft==0));
    Nsoft_1=length(find(fai_soft==1));
    Zsoft_m(timestep)=(C_soft-Nsoft_1)/(soft_num-Nsoft_1-Nsoft_0);
    Nrigid_0=length(find(fai_rigid==0));
    Nrigid_1=length(find(fai_rigid==1));
    Zrigid_m(timestep)=(C_rigid-Nrigid_1)/(rigid_num-Nrigid_1-Nrigid_0);
    Z_m(timestep)=(C_soft-Nsoft_1+C_rigid-Nrigid_1)/(particle_num-Nsoft_1-Nsoft_0-Nrigid_1-Nrigid_0);
end
% [weyy]=Axialstrain(proportion);weyy=weyy([1:10:350,350]);
% Z=Z([1:10:350,350]);Z_m=Z_m([1:10:350,350]);
% figure,plot(weyy,DDZm,'r',weyy,DRZm,'g',weyy,RRZm,'b');
% figure,plot(weyy,Z,'r',weyy,Z_m,'b');