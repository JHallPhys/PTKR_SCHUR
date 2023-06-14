% !!!!!!
% The defn of the gain and loss sets for Floquet system here appear to be
% the wrong way round. E.g gain set sorts s.t mu_1<mu_2<........
% this is because we map mu_1->1, mu_2-> 2, ..... and mu_N the largest
% positive imaginasry part to N. The index with N is the top left corner of
% new schur matrix and the N-1, N-2, etc the next prinicpal values down
% diagonal.

function [Psi_out,E_out] = REig(E_in,Psi_in,N,efn_set)

lambda=ordeig(E_in); % get order of eigs

E=-1i*log(lambda); % Calculate quasienergies23dg


if isequal(efn_set,'G') % Gain states
    
%     [Er,ind]=sort(imag(E),'descend');
    [Er,ind]=sort(real(log((lambda))))
     Er
%     [Er,ind]=sort(real(log(diag(E_in))));
    dummy(ind)=1:N;
    [Psi_out,E_out]=ordschur(Psi_in,E_in,dummy);
 
    
elseif isequal(efn_set,'S') % Stable states
    
    [Er,ind]=sort(abs(abs(diag(E_in))-1),'descend'); 
    dummy(ind)=1:N;
    [Psi_out,E_out]=ordschur(Psi_in,E_in,dummy);

        
elseif isequal(efn_set,'L') % Loss states1
    
    [Er,ind]=sort(imag(E)); 
    dummy(ind)=1:N;
    [Psi_out,E_out]=ordschur(Psi_in,E_in,dummy);
    
    
end


end