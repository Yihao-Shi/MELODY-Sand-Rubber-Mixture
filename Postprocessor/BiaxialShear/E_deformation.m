function D_strain=E_deformation(particle_num,soft_num,FRACTION)
Mass=dlmread([FRACTION 'mass_info.asc'],' ');
Jacobian=dlmread([FRACTION 'jacobian_info.asc'],' ');
rho=110000000;
E=700000;
miu=0.495;
G=E/2/(1+miu);
lamda=2*G*miu/(1-2*miu);
[softFieldNode,~]=NodeInfo(FRACTION,particle_num,soft_num);

for i=1:36
   mass(i,:)=Mass((i-1)*(sum(softFieldNode)+2)+3:i*(sum(softFieldNode)+2),1).'; 
   jacobian(i,:)=Jacobian((i-1)*(sum(softFieldNode)+2)+3:i*(sum(softFieldNode)+2),1).'; 
end
area=mass/rho;

%% deformation 
D_strain=zeros(36,1);
for timestep=1:36 
    areaPerSec=area(timestep,:).';
    jacobianPerSec=jacobian(timestep,:).';
    strainxx=dlmread([FRACTION 'GreenStrain_info/GreenXX' num2str(timestep-1) '.txt'],' ');
    strainxy=dlmread([FRACTION 'GreenStrain_info/GreenXY' num2str(timestep-1) '.txt'],' ');
    strainyy=dlmread([FRACTION 'GreenStrain_info/GreenYY' num2str(timestep-1) '.txt'],' ');
    C_xx=1+2*strainxx;C_xy=2*strainxy;C_yy=1+2*strainyy;  
    I1=C_xx+C_yy;
    psi=0.5*G*(I1-3)-G*log(jacobianPerSec)+0.5*lamda*(log(jacobianPerSec)).^2;
    D_strain(timestep)=sum(psi.*areaPerSec);
end
% for timestep=36:-1:2
%     D_strain(timestep)=D_strain(timestep)-D_strain(timestep-1);
% end
D_strain(:,1)=D_strain(:,1)-D_strain(1);
end