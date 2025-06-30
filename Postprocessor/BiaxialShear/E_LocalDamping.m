function localdamp=E_LocalDamping(particle_num,soft_num,FRACTION)
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
velxs=dlmread([FRACTION 'xvel_info.asc'],' ');
velys=dlmread([FRACTION 'yvel_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
kt=10000000000;kn=10000000000;rho=2650e5;

%% extract data
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);
for i=1:36
   master(i,:)=masters((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   len(i,:)=lens((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   velx(i,:)=velxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   vely(i,:)=velys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locx(i,:)=locxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locy(i,:)=locys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).';
end

%%
E_localdamp=zeros(36,1);localdamp=zeros(36,1);
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);
    
    velxPerSec=zeros(particle_num,max(countorNode));
    velyPerSec=zeros(particle_num,max(countorNode));
    xnormPerSec=zeros(particle_num,max(countorNode));
    ynormPerSec=zeros(particle_num,max(countorNode));
    locxPerSec=zeros(particle_num,max(countorNode));
    locyPerSec=zeros(particle_num,max(countorNode));
    lengthPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            velxPerSec(i,j)=velx(timestep,sum(countorNode(1:i-1))+j);
            velyPerSec(i,j)=vely(timestep,sum(countorNode(1:i-1))+j);
            xnormPerSec(i,j)=xnorm(timestep,sum(countorNode(1:i-1))+j);
            ynormPerSec(i,j)=ynorm(timestep,sum(countorNode(1:i-1))+j);
            locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
            locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
            lengthPerSec(i,j)=len(timestep,sum(countorNode(1:i-1))+j);
        end
    end
    
    area=zeros(particle_num,1);
    for i=1:particle_num
       x=locxPerSec(i,:);y=locyPerSec(i,:);
       for j=1:countorNode(i)-2
           area(i)=area(i)+0.5*(x(1)*y(j+1)-y(1)*x(j+1)+x(j+1)*y(j+2)-y(j+1)*x(j+2)+x(j+2)*y(1)-y(j+2)*x(1));
       end
    end

    fc=zeros(size(masterPerSec,1)*size(masterPerSec,2),19);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)~=-1
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
                if fnPerSec(i,j)~=0&&masterPerSec(i,j)<5025&&lengthPerSec(i,j)>0
                    fc(sum(countorNode(1:i-1))+j,3)=xnormPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,4)=ynormPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,5)=velxPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,6)=velyPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,7)=lengthPerSec(i,j);
                    fc(sum(countorNode(1:i-1))+j,8)=1/(rho*area(i));
                    fc(sum(countorNode(1:i-1))+j,9)=1/(rho*area(masterPerSec(i,j)+1));
                    fc(sum(countorNode(1:i-1))+j,12)=mean(velxPerSec(masterPerSec(i,j)+1,masterPerSec(masterPerSec(i,j)+1,:)==i-1));
                    fc(sum(countorNode(1:i-1))+j,13)=mean(velyPerSec(masterPerSec(i,j)+1,masterPerSec(masterPerSec(i,j)+1,:)==i-1));
                end
            end
        end
    end
    fc(fc(:,7)==0,:)=[];
    fc(isnan(fc(:,12)),:)=[];
    
    fc(:,16)=fc(:,3).*fc(:,5)+fc(:,4).*fc(:,6);
    fc(:,17)=abs(-fc(:,4).*fc(:,5)+fc(:,3).*fc(:,6));
    fc(:,18)=fc(:,3).*fc(:,12)+fc(:,4).*fc(:,13);
    fc(:,19)=abs(-fc(:,4).*fc(:,12)+fc(:,3).*fc(:,13));
    
    for i=1:size(fc,1)
        if fc(i,1)>=soft_num&fc(i,2)>=soft_num
            localdamp(timestep)=localdamp(timestep) ...
                +(abs(fc(i,16)-fc(i,18))^2*2*sqrt(1/(fc(i,8)+fc(i,9))*kn*fc(i,7)) ...
                +abs(fc(i,17)-fc(i,19))^2*2*sqrt(1/(fc(i,8)+fc(i,9))*kt*fc(i,7)));
        end
    end
end
end