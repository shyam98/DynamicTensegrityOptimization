function [M,m] = Mass(D,I_D,n_s,n_b,s,b,C_sT,C_bT,rho_s,rho_b,m_load, A_s, A_b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%  Finding the mass of each string and each bar
%  and put the mass into m_s and m_b

%Surface string mass
Number_of_string = 1;
for i =1:D:D*(n_s/2)-1
    si = sqrt(s(i,1)^2 + s(i+1,1)^2 + s(i+2,1)^2 );
    m_s(Number_of_string,1) = rho_s * A_s(1) * si ;
    Number_of_string = Number_of_string + 1 ;
end
%Internal string mass
for i =D*(n_s/2):D:D*n_s-1
    si = sqrt(s(i,1)^2 + s(i+1,1)^2 + s(i+2,1)^2 );
    m_s(Number_of_string,1) = rho_s * A_s(2) * si ;
    Number_of_string = Number_of_string + 1 ;
end


Number_of_bar = 1;
for i =1:D:D*n_b-1
    bi = sqrt(b(i,1)^2 + b(i+1,1)^2 + b(i+2,1)^2 );
    m_b(Number_of_bar,1) = rho_b * A_b * bi ;
    Number_of_bar = Number_of_bar + 1 ;
end

% Finding the mass of each node
m = 1/2 * abs(C_sT) * m_s + 1/2 * abs(C_bT) * m_b;
nnodes = size(C_sT);
m(nnodes(1),1) = m(nnodes(1),1) + m_load;

m_hat = diag(m);
M = kron(m_hat,I_D);

end

