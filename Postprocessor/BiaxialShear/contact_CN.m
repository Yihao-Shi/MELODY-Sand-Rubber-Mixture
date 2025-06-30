clear
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');

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

%%
Z=zeros(size(master,1),1);Zm=zeros(size(master,1),1);
RRZm=zeros(size(master,1),1);RRZ=zeros(size(master,1),1);
DRZm=zeros(size(master,1),1);DRZ=zeros(size(master,1),1);
DDZm=zeros(size(master,1),1);DDZ=zeros(size(master,1),1);
for timestep=1:size(master,1)
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,~]=master_process(particle_num,countorNode,master,fn,timestep);
    
    RR_num=zeros(particle_num,1);DR_num=zeros(particle_num,1);DD_num=zeros(particle_num,1);
    fc=zeros(size(masterPerSec,1)*size(masterPerSec,2),4);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)>0&&masterPerSec(i,j)<5025
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
            end
        end
    end
    fc(fc(:,1)==0&fc(:,2)==0,:)=[];

    temp=zeros(size(fc,1),7);count=1;
    [slavePar,srows]=unique(fc(:,1));srows(length(srows)+1)=size(fc,1)+1;
    for i=1:length(slavePar)
        test=sortrows(fc(srows(i):srows(i+1)-1,:),2);
        [masterPar,mrows]=unique(test(:,2));mrows(length(mrows)+1)=srows(i+1)-srows(i)+1;
        for j=1:length(masterPar)
            temp(count,1)=slavePar(i);temp(count,2)=masterPar(j);
            count=count+1;
        end
    end
    temp(count-1:end,:)=[];

    for i=1:size(temp,1)
        if temp(i,1)<soft_num&&temp(i,2)<soft_num
            DD_num(temp(i,1)+1)=DD_num(temp(i,1)+1)+1;
        elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
            DR_num(temp(i,1)+1)=DR_num(temp(i,1)+1)+1;
        else
            RR_num(temp(i,1)+1)=RR_num(temp(i,1)+1)+1;
        end
    end
    num=DD_num+DR_num+RR_num;
    
    n0=length(find(num==0));n1=length(find(num==1));
    RR0=length(find(num(soft_num+1:particle_num)==0));RR1=length(find(num(soft_num+1:particle_num)==1));
    DD0=length(find(num(1:soft_num)==0));DD1=length(find(num(1:soft_num)==1));
    RRZ(timestep)=(sum(RR_num))/(rigid_num);
    RRZm(timestep)=(sum(RR_num)-RR1)/(rigid_num-RR1-RR0);
    DRZ(timestep)=(sum(DR_num))/(particle_num);
    DRZm(timestep)=(sum(DR_num)-n1)/(particle_num-n0-n1);
    DDZ(timestep)=(sum(DD_num))/(soft_num);
    DDZm(timestep)=(sum(DD_num)-DD1)/(soft_num-DD0-DD1);
    Z(timestep)=sum(num)/particle_num;
    Zm(timestep)=(sum(num)-n1)/(particle_num-n1-n0);
end
[weyy]=Axialstrain(proportion);weyy=weyy([1:10:350,350]);
figure,plot(weyy,DDZm,'r',weyy,DRZm,'g',weyy,RRZm,'b');
figure,plot(weyy,DDZ,'r',weyy,DRZ,'g',weyy,RRZ,'b');
% figure,plot(weyy,Z,'r',weyy,Zm,'b');