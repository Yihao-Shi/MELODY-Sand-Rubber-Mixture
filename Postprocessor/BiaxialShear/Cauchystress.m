clear;
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
xpos=dlmread([FRACTION 'posx_info.asc'],' ');
ypos=dlmread([FRACTION 'posy_info.asc'],' ');
xpos=xpos(:,3:end-1);ypos=ypos(:,3:end-1);

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);

for i=1:36
   master(i,:)=masters((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).';
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

info=textread('../../Simulation_DATA/BiaxialShearTest/ball_pos_PFC/Ball&Wall_info.txt','%f');
x_min=info(1);x_max=info(2);
y_min=info(3);y_max=info(4);
Disp=dlmread([FRACTION 'displacement_info.asc'],' ');Disp=Disp([1:10:350,350],:);
lx_0=x_max-x_min-Disp(1,7)+Disp(1,9);
ly_0=y_max-y_min-Disp(1,4)+Disp(1,6);
wlx=lx_0-(Disp(:,7)-Disp(1,7))+(Disp(:,9)-Disp(1,9));
wly=ly_0-(Disp(:,4)-Disp(1,4))+(Disp(:,6)-Disp(1,6));
S=wlx.*wly;
%% 
qc=zeros(36,1);pc=zeros(36,1);
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);
    
    fxPerSec=zeros(particle_num,max(countorNode));
    xposPerSec=zeros(particle_num,max(countorNode));
    fyPerSec=zeros(particle_num,max(countorNode));
    yposPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            fxPerSec(i,j)=fx(timestep,sum(countorNode(1:i-1))+j);
            xposPerSec(i,j)=xpos(timestep,i);
            fyPerSec(i,j)=fy(timestep,sum(countorNode(1:i-1))+j);
            yposPerSec(i,j)=ypos(timestep,i);
        end
    end

    fc=zeros(size(fyPerSec,1)*size(fyPerSec,2),4);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)>0&&masterPerSec(i,j)<5025
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,3)=fxPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,4)=fyPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,5)=xposPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,1)-xposPerSec(fc(sum(countorNode(1:i-1))+j,2)+1,1);
                fc(sum(countorNode(1:i-1))+j,6)=yposPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,1)-yposPerSec(fc(sum(countorNode(1:i-1))+j,2)+1,1);
            end
        end
    end
    fc(fc(:,5)==0&fc(:,6)==0,:)=[];
    
    sigma=zeros(2,2);
    for i=1:size(fc,1)
        sigma(1,1)=sigma(1,1)+(fc(i,5)*fc(i,3))/S(timestep);
        sigma(1,2)=sigma(1,2)+(fc(i,5)*fc(i,4))/S(timestep);
        sigma(2,1)=sigma(2,1)+(fc(i,6)*fc(i,3))/S(timestep);
        sigma(2,2)=sigma(2,2)+(fc(i,6)*fc(i,4))/S(timestep);
    end
    pc(timestep)=(sigma(1,1)+sigma(2,2))/2;
    qc(timestep)=sqrt(0.5*((sigma(1,1)-pc(timestep))^2+(sigma(2,2)-pc(timestep))^2+sigma(1,2)^2+sigma(2,1)^2));
end
weyy=Axialstrain(proportion);weyy=weyy([1:10:350,350]);
plot(weyy,qc./pc);hold on;