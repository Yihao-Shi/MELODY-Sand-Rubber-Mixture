function [E_k,D_k,R_k]=E_kinetic(soft_num,particle_num,FRACTION)
energy=dlmread([FRACTION 'energy_info.asc'],' ');
energy=energy(:,3:end);
delta=zeros(36,2);energy=[energy delta];

for i=1:36
    energy(i,particle_num+1)=energy(i,particle_num+1)+sum(energy(i,1:soft_num));
    energy(i,particle_num+2)=energy(i,particle_num+2)+sum(energy(i,soft_num+1:particle_num));
    energy(i,particle_num+3)=energy(i,particle_num+3)+sum(energy(i,1:particle_num));
end
D_k=energy(:,end-2)-energy(1,end-2);
R_k=energy(:,end-1)-energy(1,end-1);
E_k=energy(:,end)-energy(1,end);
end
