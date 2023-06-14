% Function to construct the coherent states. 
% First a target directory is scanned for the coherent state.
% If one exists there then it is loaded in
% If not the function constructs the coherent state and saves to a target
% directory

function cs_out=Cs_create_component(j,norm,N,q,z,cs_in)


    
% if j==N
%     cs_in(:,:)=0;
% else
    
    
for m=-2:2
% cs_in(:,:) = cs_in(:,:)+exp(-pi*N*0.5*(abs(z).^2-z.^2)-pi*N*(q(j)-z+m).^2)*exp(1i*pi*m);
cs_in(:,:) = cs_in(:,:)+exp(-pi*N*0.5*(abs(z).^2-z.^2)-pi*N*(q(j)-z+m).^2)*exp(1i*pi*m);
end
% end
cs_out=cs_in.*norm;






end