function [D_strain,D_damp,E_contact,E_boundary,E_k]=E_work(FRACTION)

work=dlmread([FRACTION 'Work_info.asc'],' ');
kinetic=dlmread([FRACTION 'kinetic_info.asc'],' ');

%%
E_contact=work(:,4);E_contact=E_contact-E_contact(1);
D_strain=work(:,3);D_strain=D_strain-D_strain(1);
D_damp=work(:,8);D_damp=D_damp-D_damp(1);
E_boundary=sum(work(:,6:7),2);E_boundary=E_boundary-E_boundary(1);
E_k=sum(kinetic(:,3:end),2);E_k=E_k-E_k(1);    

end