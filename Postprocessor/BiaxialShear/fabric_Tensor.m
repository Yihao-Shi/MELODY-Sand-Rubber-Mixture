clear;
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
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
   len(i,:)=lens((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end
%%
a_c=zeros(size(master,1),1);a_cDD=zeros(size(master,1),1);a_cDR=zeros(size(master,1),1);a_cRR=zeros(size(master,1),1);
a_d=zeros(size(master,1),1);a_dDD=zeros(size(master,1),1);a_dDR=zeros(size(master,1),1);a_dRR=zeros(size(master,1),1);
a_n=zeros(size(master,1),1);a_nDD=zeros(size(master,1),1);a_nDR=zeros(size(master,1),1);a_nRR=zeros(size(master,1),1);
a_t=zeros(size(master,1),1);a_tDD=zeros(size(master,1),1);a_tDR=zeros(size(master,1),1);a_tRR=zeros(size(master,1),1);
bins=20;
dbin=(-90:180/bins:90)*pi/180;
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,~]=master_process(particle_num,countorNode,master,fn,timestep);

    fxPerSec=zeros(particle_num,max(countorNode));
    fyPerSec=zeros(particle_num,max(countorNode));
    xnormPerSec=zeros(particle_num,max(countorNode));
    ynormPerSec=zeros(particle_num,max(countorNode));
    xposPerSec=zeros(particle_num,max(countorNode));
    yposPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            fxPerSec(i,j)=fx(timestep,sum(countorNode(1:i-1))+j);
            fyPerSec(i,j)=fy(timestep,sum(countorNode(1:i-1))+j);
            xnormPerSec(i,j)=xnorm(timestep,sum(countorNode(1:i-1))+j);
            ynormPerSec(i,j)=ynorm(timestep,sum(countorNode(1:i-1))+j);
            xposPerSec(i,j)=xpos(timestep,i);
            yposPerSec(i,j)=ypos(timestep,i);
        end
    end
    
    fc=zeros(size(masterPerSec,1)*size(masterPerSec,2),7);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)>0&&masterPerSec(i,j)<5025
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,3)=xnormPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,4)=ynormPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,5)=fxPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,6)=fyPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,7)=abs(xposPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,1)-xposPerSec(fc(sum(countorNode(1:i-1))+j,2)+1,1));
                fc(sum(countorNode(1:i-1))+j,8)=abs(yposPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,1)-yposPerSec(fc(sum(countorNode(1:i-1))+j,2)+1,1));
            end
        end
    end
    fc(fc(:,5)==0&fc(:,6)==0,:)=[];
    
    temp=zeros(size(fc,1),11);count=1;
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
            temp(count,7)=mean(test(mrows(j):mrows(j+1)-1,7));
            temp(count,8)=mean(test(mrows(j):mrows(j+1)-1,8));
            count=count+1;
        end
    end
    temp(count-1:end,:)=[];
    temp(:,9)=temp(:,3).*temp(:,5)+temp(:,4).*temp(:,6);
    temp(:,10)=abs(-temp(:,3).*temp(:,6)+temp(:,4).*temp(:,5));
    
    RR=[];DR=[];DD=[];
    for i=1:size(temp,1)
        if temp(i,1)<soft_num&&temp(i,2)<soft_num
            DD=[DD;temp(i,3),temp(i,4),temp(i,9),temp(i,10),sqrt(temp(i,7)^2+temp(i,8)^2)];
        elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
            DR=[DR;temp(i,3),temp(i,4),temp(i,9),temp(i,10),sqrt(temp(i,7)^2+temp(i,8)^2)];
        else
            RR=[RR;temp(i,3),temp(i,4),temp(i,9),temp(i,10),sqrt(temp(i,7)^2+temp(i,8)^2)];
        end
    end
    Total=[DD;DR;RR];
    
    [acDD,acDR,acRR,ac,devia_acDD,devia_acDR,devia_acRR,devia_ac]=Fabric_ac(DD,DR,RR);
%     [adDD,adDR,adRR,ad]=Fabric_ad(DD,DR,RR,devia_acDD,devia_acDR,devia_acRR,devia_ac);
    [anDD,anDR,anRR,an]=Fabric_an(DD,DR,RR,devia_acDD,devia_acDR,devia_acRR,devia_ac);
%     [atDD,atDR,atRR,at]=Fabric_at(DD,DR,RR,devia_acDD,devia_acDR,devia_acRR,devia_ac);
    a_c(timestep)=ac;a_cDD(timestep)=acDD;a_cDR(timestep)=acDR;a_cRR(timestep)=acRR;
%     a_d(timestep)=ad;a_dDD(timestep)=adDD;a_dDR(timestep)=adDR;a_dRR(timestep)=adRR;
%     a_n(timestep)=an;a_nDD(timestep)=anDD;a_nDR(timestep)=anDR;a_nRR(timestep)=anRR;
%     a_t(timestep)=at;a_tDD(timestep)=atDD;a_tDR(timestep)=atDR;a_tRR(timestep)=atRR;
end
% A=0.5*(a_c-a_d+a_n+a_t);
weyy=Axialstrain(proportion);weyy=weyy([1:10:350,350]);
plot(weyy,a_c,'blue',weyy,a_cDD,'g',weyy,a_cDR,'r',weyy,a_cRR,'black');
% plot(weyy,A);