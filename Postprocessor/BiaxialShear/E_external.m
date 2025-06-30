function E_ex=E_external(FRACTION)
disp=dlmread([FRACTION 'displacement_info.asc'],' ');
force=dlmread([FRACTION 'force_info.asc'],' ');force=-force;
disp=[disp(:,4),disp(:,6),disp(:,7),disp(:,9)];
for i=size(disp,1):-1:1
    if i-1>0
        disp(i,:)=disp(i,:)-disp(i-1,:);
    else
        disp(i,:)=0;
    end
end
force=[force(:,4),force(:,6),force(:,7),force(:,9)];
energy=zeros(size(disp,1),5);
for i=2:size(disp,1)
    energy(i,1:4)=0.5*disp(i,:).*(force(i-1,:)+force(i,:));
end
energy(:,5)=energy(:,1)+energy(:,2)+energy(:,3)+energy(:,4);

E_ex=zeros(36,1);
for i=1:35
    E_ex(i+1)=sum(energy(1:i*10,end));
end
end