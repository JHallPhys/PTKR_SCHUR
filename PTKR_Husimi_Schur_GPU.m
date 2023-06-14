clear all
close all

parallel.gpu.enableCUDAForwardCompatibility(true) % RTX 3090 ;)


%==========================================================================
%   Quantum System and Schur Parameters
%==========================================================================

%==========================================================================
% g4=0.003
%==========================================================================
% System Parameters

N_1 =500; 
N = 2*N_1+1; % Hilbert space dimension
% return
K_class = 2.5% Classical Kicking

% cool at k=12.716,14.125
% K_class = 7.54545 % Classical Kicking
% K_class=12.6;
gamma = complex(0,0.001); % PT-strength 
T=2*pi/N; % Effective hbar
kick = K_class/T; % Quantum Kicking
hbar_eff=1/(2*pi*N);

% Schur Parameters

% eps=exp(imag(gamma)); % Tolerance parameter for stability classification
eps=1+1e-6;
% eps=1+hbar_eff/2
set_efn='G'; % Invarient Subspace: Gain ('G'), Stable ('S'), Loss ('L')
set_stability='+'; % Stability set: Gain('+'), Stable ('0'), Loss ('-')

% Check that valid set_efn and set_stability labels are used

% return
if ismember(set_efn,{'G','L'})~=1 % Sorting label check
    display('Bad Eigenfunction set label')
    return
end


if ismember(set_stability,{'+','-','0'})~=1 % Subspace label check
    display('Bad set stability label')
    return
end

%==========================================================================
%   Matrix Construction and Schur decompesition
%==========================================================================

U=zeros(N,N); % Initialise Flouqet matrix
[U,time1]=UMatrix(U,N,N_1,K_class,T,gamma); % Construct Flouqet matrix
[psi,En] = schur(U); % psi are the Schur eigenfns and En matrix of eigs

% return
%==========================================================================
%   Project onto Subspace of stability
%==========================================================================

[psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values

Es=diag(Es);
figure(1)
hold on
for j = 1:N
plot(j,abs(Es(j)),'k.','Markersize',1)
end
[psi_2,n_efn]=Psi_lifetime(psiS,Es,eps,set_stability);
% 
% for j=1:n_efn
%     
%     G=real(log(psi_2(:,j)'*U*psi_2(:,j)));
%     
%     figure(2)
%     hold on 
%     plot(j,G,'b.','markersize',1)
% 
% end
% 
% 
% return



n_efn
% % return
n_efn/N
% return


% ylim([-0.15 0.15])
% ylabel('Im(\epsilon)')
% xlabel('Re(\epsilon)/ \pi')
% return 
%==========================================================================
%  Arrays independent of splitying
%==========================================================================

% Discrete Phase Space Grid

q = linspace(0,1,N); % q interval
p = linspace(-1/2,1/2,N); % p interval
[qmesh,pmesh]=meshgrid(q,p); 
z = (qmesh+1i*(pmesh)); 
Hus = zeros(N,N,n_efn); % Average Husimi function array
cs=zeros(N,N);
norm_cs= (2/N)^0.25; % Normalisation constant for the coherent state
%==========================================================================
%   Create the split interval
%==========================================================================
n_split=4;
split=round(linspace(1,n_efn,n_split));
% return

% Gpu arrays independent of split;

psi_2_gpu=gpuArray(psi_2);
cs_gpu=gpuArray(cs);
q_gpu=gpuArray(q);
z_gpu=gpuArray(z);


%==========================================================================
%   Begin split Husimi
%==========================================================================

for sn=2:length(split)
    
nstart=split(sn-1);
if sn>2 
    nstart=nstart+1; 
end
nend=split(sn);

% GPU Array Decleration
ds=abs(nend-nstart)+1;
Hus_split_gpu=gpuArray(zeros(N,N,ds));
psi_split_gpu=gpuArray(zeros(1,1,ds));

tic
for j = 1:N-1
    
   
    disp([num2str(j),' out of ',num2str(N),' for ',num2str(sn-1),' out of ',num2str(length(split)-1)]) % keep track
    cs_gpu=Cs_create_component(j,norm_cs,N,q_gpu,z_gpu,cs_gpu); 
    psi_split_gpu(1,1,:)=psi_2_gpu(j,nstart:nend);
    Hus_split_gpu=Hus_split_gpu+conj(cs_gpu).*psi_split_gpu;
    cs_gpu(:,:)=0;
    
end
toc

Hus(:,:,nstart:nend)=gather(Hus_split_gpu);
clear Hus_split_gpu psi_split_gpu

end
tic
% 
% n_efn=1

Hus_av=zeros(N,N);
for t=1:n_efn
        disp([num2str(t),' out of ',num2str(n_efn)]) % keep track 
        Hus_av=Hus_av+abs(Hus(:,:,t)).^2;
%         figure(1)
%         clf
%         imagesc(q,p,Hus_av)
%         colorbar
%         caxis([0 1])
%         colormap(viridis)
%         set(gca,'YDir','normal')
%         pause(0.1)
end
time2=toc;



viridis=viridis();

figure
clf
imagesc(q,p,Hus_av)
colorbar
caxis([0 1])
colormap(viridis)
set(gca,'YDir','normal')




% figure(5)
% clf
% imagesc(q,p,Hus_av)
% colormap(parula)
% set(gca,'YDir','normal')
% caxis([0 1])


if max(max(Hus_av))>1
    'Hus_av has elements that are > 1!!!'
end


if sum(sum(isnan(Hus_av)))>0
    'Hus_av has elements that are NaN!!!'
end
% return
% while 1 
% 
% % Tell the user about the system and ask for input
%     
% user_msg_0 = strcat('Max Efns: ',num2str(n_efn));    
% display(user_msg_0);
% user_msg_1=' Select the number of eigenfunctions to look at ';
% nsplit=input(user_msg_1);
% 
% 
% Hus_av=zeros(N,N);
% for t=1:nsplit
%         disp([num2str(t),' out of ',num2str(n_efn)]) % keep track 
%         Hus_av=Hus_av+abs(Hus(:,:,t)).^2;
% end
% time2=toc;
% 
% 
% 
% viridis=viridis();
% 
% figure
% clf
% imagesc(q,p,Hus_av)
% colorbar
% % caxis([0 0.2])
% colormap(viridis)
% set(gca,'YDir','normal')
% 
% 
% user_exit=input(' Repeat with new nsplit? [y/n]:  ','s');
% 
% if isequal(user_exit,'n')
% 
%     break
%     
% end






% end












