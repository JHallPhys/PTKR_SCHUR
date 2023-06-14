% Construct the Floquet operator of size(NxN) for the kicked rotor.

function [U_out,time_out]=UMatrix(U_in,N,N_1,K_class,T,gamma)

K_s = K_class/T; % Quantum Kicking

j_i=linspace(-N_1,N_1,N);
j_i=transpose(j_i);
l1=linspace(1,N,N);
tic
for k = 1:N
           
            b=exp(-1i*T*0.5*j_i.^2+imag(gamma)*j_i+2*pi*1i*j_i.*(k-l1)/N);   
            U_in(k,:) =  (exp((-1i*K_s*cos(2*pi*k/N)))*sum(b))/N;  
                
end

U_out=U_in;
time_out=toc;


end