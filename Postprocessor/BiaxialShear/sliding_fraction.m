clear;
proportion='10';
miu=0.5;
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
   locx(i,:)=locxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locy(i,:)=locys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).';
end

%%
DD=zeros(36,1);DR=zeros(36,1);RR=zeros(36,1);sf=zeros(36,1);wf=zeros(36,1);weak=cell(1,1);strong=cell(1,1);
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);

    fxPerSec=zeros(particle_num,max(countorNode));
    fyPerSec=zeros(particle_num,max(countorNode));
    xposPerSec=zeros(particle_num,max(countorNode));
    yposPerSec=zeros(particle_num,max(countorNode));
    locxPerSec=zeros(particle_num,max(countorNode));
    locyPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            fxPerSec(i,j)=fx(timestep,sum(countorNode(1:i-1))+j);
            fyPerSec(i,j)=fy(timestep,sum(countorNode(1:i-1))+j);
            xposPerSec(i,j)=xpos(timestep,i);
            yposPerSec(i,j)=ypos(timestep,i);
            locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
            locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
            xnormPerSec(i,j)=xnorm(timestep,sum(countorNode(1:i-1))+j);
            ynormPerSec(i,j)=ynorm(timestep,sum(countorNode(1:i-1))+j);
        end
    end
    
    area=zeros(particle_num,1);
    for i=1:particle_num
       x=locxPerSec(i,:);y=locyPerSec(i,:);
       for j=1:countorNode(i)-2
           area(i)=area(i)+0.5*(x(1)*y(j+1)-y(1)*x(j+1)+x(j+1)*y(j+2)-y(j+1)*x(j+2)+x(j+2)*y(1)-y(j+2)*x(1));
       end
    end

    DDcount=0;DRcount=0;RRcount=0;
    fc=zeros(size(fxPerSec,1)*size(fxPerSec,2),12);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)~=-1
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,3)=xnormPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,4)=ynormPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,5)=fxPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,6)=fyPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,7)=locxPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,j)-xposPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,1);
                fc(sum(countorNode(1:i-1))+j,8)=locyPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,j)-yposPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,1);
                dist=sqrt(fc(sum(countorNode(1:i-1))+j,7)^2+fc(sum(countorNode(1:i-1))+j,8)^2);
                fc(sum(countorNode(1:i-1))+j,9)=fc(sum(countorNode(1:i-1))+j,5)*fc(sum(countorNode(1:i-1))+j,7)/area(i);
                fc(sum(countorNode(1:i-1))+j,10)=fc(sum(countorNode(1:i-1))+j,5)*fc(sum(countorNode(1:i-1))+j,8)/area(i);
                fc(sum(countorNode(1:i-1))+j,11)=fc(sum(countorNode(1:i-1))+j,6)*fc(sum(countorNode(1:i-1))+j,8)/area(i);
                fc(sum(countorNode(1:i-1))+j,12)=fc(sum(countorNode(1:i-1))+j,6)*fc(sum(countorNode(1:i-1))+j,7)/area(i);
                fc(sum(countorNode(1:i-1))+j,13)=locxPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,14)=locyPerSec(i,j);
            end
        end
    end
    fc(fc(:,1)==0&fc(:,2)==0|fc(:,2)==5025|fc(:,2)==5026|fc(:,2)==5027|fc(:,2)==5028,:)=[];
    fc(:,15)=fc(:,3).*fc(:,5)+fc(:,4).*fc(:,6);
    fc(:,16)=-fc(:,4).*fc(:,5)+fc(:,3).*fc(:,6);
    
    temp=zeros(size(fc,1),12);count=1;
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
            temp(count,7)=sum(test(mrows(j):mrows(j+1)-1,9));
            temp(count,8)=(sum(test(mrows(j):mrows(j+1)-1,10))+sum(test(mrows(j):mrows(j+1)-1,11)))/2;
            temp(count,9)=sum(test(mrows(j):mrows(j+1)-1,12));
            count=count+1;
        end
    end
    temp(count-1:end,:)=[];
    temp(:,10)=temp(:,3).*temp(:,5)+temp(:,4).*temp(:,6);
    temp(:,11)=abs(-temp(:,3).*temp(:,6)+temp(:,4).*temp(:,5));
    temp(temp(:,1)>5023|temp(:,2)>5023,:)=[];
    
    [partNum,prows]=unique(temp(:,1));prows(length(prows)+1)=size(temp,1)+1;
    partList=zeros(size(partNum,1),7);
    for i=1:length(partNum)
        partList(i,1)=partNum(i);
        partList(i,2)=sum(temp(prows(i):prows(i+1)-1,3));
        partList(i,3)=sum(temp(prows(i):prows(i+1)-1,4));
        partList(i,4)=sum(temp(prows(i):prows(i+1)-1,5));
        partList(i,5)=0.5*(partList(i,2)+partList(i,4)+sqrt((partList(i,2)-partList(i,4))^2+4*partList(i,3)^2));
        partList(i,6)=0.5*(partList(i,2)+partList(i,4)-sqrt((partList(i,2)-partList(i,4))^2+4*partList(i,3)^2));
        partList(i,7)=0.5*atan((2*partList(i,3))/(partList(i,2)-partList(i,4)))+pi/2;
    end
    mean_force=mean(partList(:,5));
    
    i=1;
    while i<=length(partList)
        if i~=partList(i,1)+1
            partList=[partList(1:i-1,:);0,0,0,0,0,0,0;partList(i:length(partList),:)];
        end
        i=i+1;
    end
    
    for i=1:length(temp)
        masterid=temp(i,1)+1;slaveid=temp(i,2)+1;
        if partList(masterid,5)>mean_force&partList(slaveid,5)>mean_force
            temp(i,12)=1;
        end
    end

%     RR_force=[];DR_force=[];DD_force=[];
%     for i=1:size(temp,1)
%         if temp(i,1)<soft_num&&temp(i,2)<soft_num
%             DD_force=[DD_force;temp(i,8)/(miu*temp(i,7))];
%         elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
%             DR_force=[DR_force;temp(i,8)/(miu*temp(i,7))];
%         else
%             RR_force=[RR_force;temp(i,8)/(miu*temp(i,7))];
%         end
%     end

    s=[];w=[];meanf=mean(temp(:,7));
    for i=1:size(temp,1)
        ori=atan((ypos(timestep,temp(i,1)+1)-ypos(timestep,temp(i,2)+1))/(xpos(timestep,temp(i,1)+1)-xpos(timestep,temp(i,2)+1)));
        if temp(i,12)>0
            s=[s;temp(i,11)/(miu*temp(i,10)),ori];
        else
            w=[w;temp(i,11)/(miu*temp(i,10)),ori];
        end
    end
    sf(timestep)=length(find(s(:,1)>0.999))/(length(s)+length(w));
    wf(timestep)=length(find(w(:,1)>0.999))/(length(w)+length(s));
    strong{timestep,1}=s;weak{timestep,1}=w;
    
%     DD(timestep)=length(find(DD_force>0.999))/length(DD_force);
%     DR(timestep)=length(find(DR_force>0.999))/length(DR_force);
%     RR(timestep)=length(find(RR_force>0.999))/length(RR_force);
end
weyy=Axialstrain(proportion);weyy=weyy([1:10:350,350]);