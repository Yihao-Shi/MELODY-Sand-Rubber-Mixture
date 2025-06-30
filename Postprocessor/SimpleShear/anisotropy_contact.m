clear
proportion='40';
timestep=241;
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
xnorm=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorm=dlmread([FRACTION 'ynorm_info.asc'],' ');
fx=dlmread([FRACTION 'fx_info.asc'],' ');
fy=dlmread([FRACTION 'fy_info.asc'],' ');
master=dlmread([FRACTION 'master_info.asc'],' ');
len=dlmread([FRACTION 'length_info.asc'],' ');
id=[];
for i=3:size(master,2)-1
    if mod(i-2,64)==0
        continue;
    else
        id=[id,i];
    end
end
master=master(:,id);
len=len(:,id);
xnorm=xnorm(:,id);
ynorm=ynorm(:,id);
fx=fx(:,id);
fy=fy(:,id);

%% extract data
nodes=63;
particle_num=length(id)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

%%
bins=20;
dbin=(-90:180/bins:90)*pi/180;
force_direction=zeros(2*bins,2);f_n=zeros(2*bins,2);

fn=fx.*xnorm+fy.*ynorm;
[masterPerSec,~]=master_process(particle_num,nodes,master,fn,timestep);
[contact_info,~,~,~]=ContactPairs(masterPerSec,soft_num);

fxPerSec=zeros(particle_num,nodes);
fyPerSec=zeros(particle_num,nodes);
xnormPerSec=zeros(particle_num,nodes);
ynormPerSec=zeros(particle_num,nodes);
for i=1:particle_num
    for j=1:nodes
        fxPerSec(i,j)=fx(timestep,(i-1)*nodes+j);
        fyPerSec(i,j)=fy(timestep,(i-1)*nodes+j);
        xnormPerSec(i,j)=xnorm(timestep,(i-1)*nodes+j);
        ynormPerSec(i,j)=ynorm(timestep,(i-1)*nodes+j);
    end
end

for i=1:size(masterPerSec,1)
   for j=1:size(masterPerSec,2)
       end1=min(i-1,masterPerSec(i,j));
       end2=max(i-1,masterPerSec(i,j));
       for k=1:size(contact_info,1)
           if end1==contact_info(k,1)&end2==contact_info(k,2)
               if i-1<masterPerSec(i,j)
                    contact_info(k,3)=contact_info(k,3)+fxPerSec(i,j);
                    contact_info(k,4)=contact_info(k,4)+fyPerSec(i,j);
                    contact_info(k,5)=contact_info(k,5)+xnormPerSec(i,j);
                    contact_info(k,6)=contact_info(k,6)+ynormPerSec(i,j);
                    normalize=sqrt(contact_info(k,5)^2+contact_info(k,6)^2);
                    contact_info(k,5)=contact_info(k,5)/normalize;
                    contact_info(k,6)=contact_info(k,6)/normalize;
                else
                    contact_info(k,7)=contact_info(k,7)+fxPerSec(i,j);
                    contact_info(k,8)=contact_info(k,8)+fyPerSec(i,j);
                    contact_info(k,9)=contact_info(k,9)+xnormPerSec(i,j);
                    contact_info(k,10)=contact_info(k,10)+ynormPerSec(i,j);
                    normalize=sqrt(contact_info(k,9)^2+contact_info(k,10)^2);
                    contact_info(k,9)=contact_info(k,9)/normalize;
                    contact_info(k,10)=contact_info(k,10)/normalize;
                end
               break;
           end
       end
   end
end
contact_info(:,11)=contact_info(:,3).*contact_info(:,5)+contact_info(:,4).*contact_info(:,6);
contact_info(:,12)=contact_info(:,7).*contact_info(:,9)+contact_info(:,8).*contact_info(:,10);
contact_info(:,13)=abs(-contact_info(:,3).*contact_info(:,6)+contact_info(:,4).*contact_info(:,5));
contact_info(:,14)=abs(-contact_info(:,7).*contact_info(:,10)+contact_info(:,8).*contact_info(:,9));
thetaPerSec=[atan(contact_info(:,6)./contact_info(:,5));atan(contact_info(:,10)./contact_info(:,9))];
fnormPerSec=[contact_info(:,11);contact_info(:,12)];

for j=1:2*bins
    if j>bins
        force_direction(j,1)=force_direction(j-bins,1)+pi;
        f_n(j,1)=f_n(j-bins,1)+pi;
        force_direction(j,2)=force_direction(j-bins,2);
        f_n(j,2)=f_n(j-bins,2);
    else
        force_direction(j,1)=(dbin(j)+dbin(j+1))/2;
        f_n(j,1)=(dbin(j)+dbin(j+1))/2;

        row=find(thetaPerSec>dbin(j)&thetaPerSec<=dbin(j+1));
        force_direction(j,2)=length(row);
        norm=0;tan=0;
        for n=1:length(row)
            norm=norm+fnormPerSec(row(n));
        end
        f_n(j,2)=norm/length(row);
    end
end
normDirection=[contact_info(:,5),contact_info(:,6);contact_info(:,9),contact_info(:,10)];
% for i=1:size(force_direction,1)-1
%     normDirection(sum(force_direction(1:i-1,2))+1:sum(force_direction(1:i,2)),1)=cos(force_direction(i,1));
%     normDirection(sum(force_direction(1:i-1,2))+1:sum(force_direction(1:i,2)),2)=sin(force_direction(i,1));
% end
force_direction(:,2)=force_direction(:,2)/(2*pi/bins*sum(force_direction(1:end/2,2)));
%% Fitting
delta=[1,0;0,1];fai=zeros(2,2);gamma_n=zeros(2,2);
for i=1:length(normDirection)
    fai=fai+kron(normDirection(i,:),normDirection(i,:).');
end
fai=fai./length(normDirection);
para=0.5*(fai(1,1)+fai(2,2));
devia_ac=4*(fai-para*delta);
for i=1:length(normDirection)
    temp_nn=kron(normDirection(i,:),normDirection(i,:).');
    gamma_n=gamma_n+(fnormPerSec(i)*temp_nn)/(1+sum(sum(devia_ac.*temp_nn)));
end
gamma_n=gamma_n./length(normDirection);
para_n=0.5*(gamma_n(1,1)+gamma_n(2,2));
devia_an=4*(gamma_n-para_n*delta)/(gamma_n(1,1)+gamma_n(2,2));
fittingBin=(0:pi/40:2*pi).';E_xx=zeros(length(fittingBin),1);F_xx=zeros(length(fittingBin),1);
for i=1:length(fittingBin)
    n=[cos(fittingBin(i)),sin(fittingBin(i))];
    nxn=kron(n,n.');
    E_xx(i)=(1+sum(sum(devia_ac.*nxn)))/(2*pi);
    F_xx(i)=(gamma_n(1,1)+gamma_n(2,2))*(1+sum(sum(devia_an.*nxn)));
end
%% plot
figure,polarplot(force_direction(:,1),force_direction(:,2),'o',fittingBin,E_xx,'-');
figure,polarplot(f_n(:,1),f_n(:,2),'o',fittingBin,F_xx,'-');