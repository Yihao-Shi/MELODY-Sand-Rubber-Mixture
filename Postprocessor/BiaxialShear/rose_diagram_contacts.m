clear;
proportion='30';
timestep=36;
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');
xpos=dlmread([FRACTION 'posx_info.asc'],' ');
ypos=dlmread([FRACTION 'posy_info.asc'],' ');
xpos=xpos(timestep,3:end-1).';ypos=ypos(timestep,3:end-1).';

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);

for i=1:36
   master(i,:)=masters((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   len(i,:)=lens((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locx(i,:)=locxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locy(i,:)=locys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

%%
bins=20;
dbin=(-90:180/bins:90)*pi/180;force_direction=zeros(2*bins+1,2);f_n=zeros(2*bins+1,2);
fn=fx.*xnorm+fy.*ynorm;
[masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);

fxPerSec=zeros(particle_num,max(countorNode));
fyPerSec=zeros(particle_num,max(countorNode));
xnormPerSec=zeros(particle_num,max(countorNode));
ynormPerSec=zeros(particle_num,max(countorNode));
locxPerSec=zeros(particle_num,max(countorNode));
locyPerSec=zeros(particle_num,max(countorNode));
for i=1:particle_num
    for j=1:countorNode(i)
        fxPerSec(i,j)=fx(timestep,sum(countorNode(1:i-1))+j);
        fyPerSec(i,j)=fy(timestep,sum(countorNode(1:i-1))+j);
        xnormPerSec(i,j)=xnorm(timestep,sum(countorNode(1:i-1))+j);
        ynormPerSec(i,j)=ynorm(timestep,sum(countorNode(1:i-1))+j);
        locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
        locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
    end
end

fc=zeros(size(masterPerSec,1)*size(masterPerSec,2),6);
for i=1:particle_num
    for j=1:countorNode(i)
        if masterPerSec(i,j)>0&&masterPerSec(i,j)<5025
            fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
            fc(sum(countorNode(1:i-1))+j,3)=xnormPerSec(i,j);
            fc(sum(countorNode(1:i-1))+j,4)=ynormPerSec(i,j);
            fc(sum(countorNode(1:i-1))+j,5)=fxPerSec(i,j);
            fc(sum(countorNode(1:i-1))+j,6)=fyPerSec(i,j);
        end
    end
end
fc(fc(:,5)==0&fc(:,6)==0,:)=[];

temp=zeros(size(fc,1),7);count=1;
[slavePar,srows]=unique(fc(:,1));srows(length(srows)+1)=size(fc,1)+1;
for i=1:length(slavePar)
    test=sortrows(fc(srows(i):srows(i+1)-1,:),2);
    [masterPar,mrows]=unique(test(:,2));mrows(length(mrows)+1)=srows(i+1)-srows(i)+1;
    for j=1:length(masterPar)
        temp(count,1)=slavePar(i);temp(count,2)=masterPar(j);
        temp(count,3)=sum(test(mrows(j):mrows(j+1)-1,3));
        temp(count,4)=sum(test(mrows(j):mrows(j+1)-1,4));
        normalize=sqrt(temp(count,3)^2+temp(count,4)^2);
        temp(count,3)=temp(count,3)/normalize;
        temp(count,4)=temp(count,4)/normalize;
        temp(count,5)=sum(test(mrows(j):mrows(j+1)-1,5));
        temp(count,6)=sum(test(mrows(j):mrows(j+1)-1,6));
        count=count+1;
    end
end
temp(count-1:end,:)=[];
temp(:,7)=abs(temp(:,3).*temp(:,5)+temp(:,4).*temp(:,6));

RR=[];DR=[];DD=[];
for i=1:size(temp,1)
    if temp(i,1)<soft_num&&temp(i,2)<soft_num
        DD=[DD;temp(i,3),temp(i,4),temp(i,7),atan(temp(i,4)/temp(i,3))];
    elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
        DR=[DR;temp(i,3),temp(i,4),temp(i,7),atan(temp(i,4)/temp(i,3))];
    else
        RR=[RR;temp(i,3),temp(i,4),temp(i,7),atan(temp(i,4)/temp(i,3))];
    end
end
Total=[DD;DR;RR];
if proportion=='0'
    DD=[0,0;0,0];DR=[0,0;0,0];
end

for j=1:2*bins+1
    if j>bins
        force_direction(j,1)=force_direction(j-bins,1)+pi;
        f_n(j,1)=f_n(j-bins,1)+pi;
        force_direction(j,2:5)=force_direction(j-bins,2:5);
        f_n(j,2:5)=f_n(j-bins,2:5);
    else
        force_direction(j,1)=(dbin(j)+dbin(j+1))/2;
        f_n(j,1)=(dbin(j)+dbin(j+1))/2;

        force_direction(j,2)=sum(sum(Total(:,4)>dbin(j)&Total(:,4)<=dbin(j+1)));
        force_direction(j,3)=sum(sum(DD(:,4)>dbin(j)&DD(:,4)<=dbin(j+1)));
        force_direction(j,4)=sum(sum(DR(:,4)>dbin(j)&DR(:,4)<=dbin(j+1)));
        force_direction(j,5)=sum(sum(RR(:,4)>dbin(j)&RR(:,4)<=dbin(j+1)));

        f_n(j,2)=sum(Total(Total(:,4)>dbin(j)&Total(:,4)<=dbin(j+1),3))/force_direction(j,2);
        f_n(j,3)=sum(DD(DD(:,4)>dbin(j)&DD(:,4)<=dbin(j+1),3))/force_direction(j,3);
        f_n(j,4)=sum(DR(DR(:,4)>dbin(j)&DR(:,4)<=dbin(j+1),3))/force_direction(j,4);
        f_n(j,5)=sum(RR(RR(:,4)>dbin(j)&RR(:,4)<=dbin(j+1),3))/force_direction(j,5);
    end
end
for c=2:5
    force_direction(:,c)=force_direction(:,c)/(2*pi/bins*sum(force_direction(1:(end-1)/2,c)));
end
%% Fitting
[~,~,~,~,devia_acDD,devia_acDR,devia_acRR,devia_ac]=Fabric_ac(DD,DR,RR);
[~,~,~,~,gamma_n,gammaDD_n,gammaDR_n,gammaRR_n,devia_anDD,devia_anDR,devia_anRR,devia_an]=Fabric_an(DD,DR,RR,devia_acDD,devia_acDR,devia_acRR,devia_ac);
fittingBin=(0:pi/40:2*pi).';E_xx=zeros(length(fittingBin),5);F_xx=zeros(length(fittingBin),5);
for i=1:length(fittingBin)
    n=[cos(fittingBin(i)),sin(fittingBin(i))];
    nxn=kron(n,n.');
    E_xx(i,1)=fittingBin(i);F_xx(i,1)=fittingBin(i);
    E_xx(i,2)=(1+sum(sum(devia_ac.*nxn)))/(2*pi);
    E_xx(i,3)=(1+sum(sum(devia_acDD.*nxn)))/(2*pi);
    E_xx(i,4)=(1+sum(sum(devia_acDR.*nxn)))/(2*pi);
    E_xx(i,5)=(1+sum(sum(devia_acRR.*nxn)))/(2*pi);
    F_xx(i,2)=(gamma_n(1,1)+gamma_n(2,2))*(1+sum(sum(devia_an.*nxn)));
    F_xx(i,3)=(gammaDD_n(1,1)+gammaDD_n(2,2))*(1+sum(sum(devia_anDD.*nxn)));
    F_xx(i,4)=(gammaDR_n(1,1)+gammaDR_n(2,2))*(1+sum(sum(devia_anDR.*nxn)));
    F_xx(i,5)=(gammaRR_n(1,1)+gammaRR_n(2,2))*(1+sum(sum(devia_anRR.*nxn)));
end
%% plot
for c=2:5
    figure,polarplot(force_direction(:,1),force_direction(:,c),'o',E_xx(:,1),E_xx(:,c),'-');
    figure,polarplot(f_n(:,1),f_n(:,c),'o',F_xx(:,1),F_xx(:,c),'-');
end