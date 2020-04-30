function [N,Cb,Cs,nnodes,n_s,n_b, zl_i] = Lander_3D(q,p,r,L,cyl, C_2, z_position)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% TENSEGRITY GEOMETRY TORUS
% p: Number of side 
% p: Number of level

%L = 3; 

%cyl = 'SP';
% 'RCC': right circular cylinder
% 'REC': right elliptic cylinder
% 'POR': paraboloid of revolution
% 'SP': sphere

% RCC 
% REC
ar = 0.9;
br = 0.5;
% POR
Ld = L+0.01*L; % > L
a = 0.1;
% SP
Rd = L/2 + 10^(-2);%0.5*(L+0.01*L); % > L/2

figure(1)
%axis([-0.75 0.75 -0.75 0.75 -0.25 1.25])

% VIEW
az = 0;
el = 90; % 0:y, 90:z


% CREATING N MATRIX
N = zeros(3,2*p*q+p);
%
c = 1; % counter



for i =0:q
    for k = 0:p-1
        for l = 1:2
            dontp = 0;
            %Linear
            %zl_p = (2*i + (l-1))
            %Second Order
            zl_p = (C_2 + 2*C_2*q-1 + (1-2*C_2*q-2*C_2)*(2*i+l) + C_2*(2*i+l)^2);
            
            zl = zl_p*L/(2*q);
           
            zl_i(i+1) = zl;
            zlv = [0;0;zl];
            tetl = (2*k + l)*pi/p;
            
            switch cyl
                case 'RCC' % 'RCC': right circular cylinder
                    ar = r;
                    br = r;
                    center_node = [0;0;z_position*L];
                case 'REC' % 'REC': right elliptic cylinder
                    
                case 'POR' % 'POR': paraboloid of revolution
                    %ar = (4*a*(Ld - zl))^(1/2);
                    %br = (4*a*(Ld - zl))^(1/2);
                    ar = (4*a*(Ld - (-zl + L)))^(1/2);
                    br = (4*a*(Ld - (-zl + L)))^(1/2);
                    center_node = [0;0;z_position*L];
                case 'SP' % 'SP': sphere
                    ar = (Rd^2 - (zl-L/2)^2)^(1/2);
                    br = (Rd^2 - (zl-L/2)^2)^(1/2);
                    center_node = [0;0;z_position*r];
                otherwise
                    disp('Wrong shape input')
                    disp('Exiting')
                    return
            end
            rl = [ar*cos(tetl); br*sin(tetl); 0];        
            if i == q
                if l == 2
                    dontp = 1;
                end
            end
            if dontp == 0
                N(1:end,c) = zlv + rl;
                c = c + 1;
            end
        end
    end
end
N(3,:) = N(3,:) - L/2;

%
%
figure(1)
hold on
l = 1;
for i = 1:2*p*q+p
    if i > 2*p*q
        plot3(N(1,i),N(2,i),N(3,i),'ok','LineWidth',1,'MarkerSize',4)
    else
        if l == 1
            l = 2;
            plot3(N(1,i),N(2,i),N(3,i),'ok','LineWidth',0.5,'MarkerSize',4)
        else
            l = 1;
            plot3(N(1,i),N(2,i),N(3,i),'sk','LineWidth',0.5,'MarkerSize',4)
        end
    end
end

xlabel('x'); ylabel('y'); zlabel('z');
%grid on
view(az, el);


%
% CONNECTIVITY BARS
E1 = [1 0; 0 0];
E2 = [0 0; 0 1];
PHI = zeros(2*p);
PHIb = zeros(2*p,3*p);
PHIt = zeros(p,2*p);
IIb = zeros(2*p,3*p);

e1 = [1 0]'; e2 = [0 1]'; z0 = [0 0]';

E3b = [zeros(2,2), e2];
E4b = [eye(2,2), e1];
E2b = [E2, z0];
E1b = [E1, z0];
E2t = e2';
E1t = e1';

for i=1:p-1
    i1 = 2*(i-1)+1;
    i2 = 2*i+1;
    PHI(i1:i1+1, i2:i2+1) = E2;
    PHI(i2:i2+1, i1:i1+1) = E1;
end
PHI(1:2,end-1:end) = E1;
PHI(end-1:end,1:2) = E2;

for i=1:p-1
    i1x = 2*(i-1)+1;
    i1y = 3*i+1;
    i2x = 2*i+1;
    i2y = 3*(i-1)+1;    
    PHIb(i1x:i1x+1, i1y:i1y+2) = E2b;
    PHIb(i2x:i2x+1, i2y:i2y+2) = E1b;
end
PHIb(1:2,end-2:end) = E1b;
PHIb(end-1:end,1:3) = E2b;

for i=1:p
    i1 = 2*(i-1)+1;
    i2 = 3*(i-1)+1;
    IIb(i1:i1+2-1, i2:i2+3-1) = E4b;

    ii1 = (i-1)+1;
    PHIt(ii1,i1:i1+1) = E2t;
    
    if i ~=p
        i2 = 3*i+1;
        ii1 = ii1+1;
        IIb(i1:i1+2-1, i2:i2+3-1) = -E3b;
        PHIt(ii1,i1:i1+1) = E1t;
    end
end
IIb(end-2+1:end, 1:3) = -E3b;
PHIt(1,end-2+1:end) = E1t;

CB = zeros((2*q+1)*p);
CB(1:2*p,1:3*p) = -IIb;
CB(2*p+1:2*p+1+2*p-1,1:3*p) = PHIb;
for i = 1:q-1
    i1x = 2*p + 2*p*(i-1)+1;
    i1y = 3*p + 2*p*(i-1)+1;
    i2 = 2*p + 2*p*i + 1;
    
    CB(i1x:i1x+2*p-1, i1y:i1y+2*p-1) = -eye(2*p);
    if i ~= q-1
        CB(i2:i2+2*p-1, i1y:i1y+2*p-1) = PHI;
    end
end

CB(end-p+1:end,end-2*p+1:end) = PHIt;
%
CB = CB';
%
B = N*CB';

CBm1 = zeros(size(CB,1));% CB-1
CBp1 = zeros(size(CB,2));
for i = 1:size(CB,1)
    for j = 1:size(CB,2)
        if CB(i,j) == -1
            CBm1(i,j) = -CB(i,j);
        elseif CB(i,j) == 1
            CBp1(i,j) = CB(i,j);
        end
    end
end

Bm1 = N*CBm1';
Bp1 = N*CBp1';



%PLOTTING BARS
l=1;
for i = 1:size(B,2)
    if l == 1
        l = 2;
        plot3([Bm1(1,i),Bp1(1,i)], ...
            [Bm1(2,i),Bp1(2,i)], ...
            [Bm1(3,i),Bp1(3,i)],...
            'b','LineWidth',1.75)
    else
        l = 1;
        plot3([Bm1(1,i),Bp1(1,i)], ...
            [Bm1(2,i),Bp1(2,i)], ...
            [Bm1(3,i),Bp1(3,i)],...
            'b','LineWidth',1.75)
    end

end


%
% % STRING CONNECTIVITY

TET1 = [z0 z0 e2 -e2 -e2 z0 z0 -e2];
TET2 = [e2-e1 -e2 z0 e1 e2 -e1 z0 z0];
TET3 = [zeros(2,6), -e1 e2];
TET4 = [z0 e1 -e1 z0 z0 e1 e1 z0];
PSI1 = zeros(2*p,8*p);
PSI2 = zeros(2*p,8*p);

PSI1b = zeros(2*p,8*p);
PSI2b = zeros(2*p,8*p);
PSI1t = zeros(2*p,6*p);
PSI2t = zeros(p,6*p);

TET1b = [TET1(1:end,1:3), TET1(1:end,5:8), -e1];
TET2b = [TET2(1:end,1:3), TET2(1:end,5:8), e1];
TET3b = [TET3(1:end,1:3), TET3(1:end,5:8), z0];
TET4b = [TET4(1:end,1:3), TET4(1:end,5:8), z0];

TET1t = [TET1(1:end,1), TET1(1:end,3:7)];
TET2t = [TET2(1:end,1), TET2(1:end,3:7)];
TET3t = [1, 0]*[TET3(1:end,1), TET3(1:end,3:7)];
TET4t = [1, 0]*[TET4(1:end,1), TET4(1:end,3:7)];
for i=1:p
    i1 = 2*(i-1)+1;
    i2 = 8*(i-1)+1;
    i2t = 6*(i-1)+1;
    i1t = (i-1)+1;
    PSI1(i1:i1+2-1, i2:i2+8-1) = TET2;
    PSI2(i1:i1+2-1, i2:i2+8-1) = TET4;
    PSI1b(i1:i1+2-1, i2:i2+8-1) = TET2b;
    PSI2b(i1:i1+2-1, i2:i2+8-1) = TET4b;
    PSI1t(i1:i1+2-1, i2t:i2t+6-1) = TET2t;
    PSI2t(i1t, i2t:i2t+6-1) = TET4t;    
    if i ~=p
        i2 = 8*i+1;
        i2t = 6*i+1;
        PSI1(i1:i1+2-1, i2:i2+8-1) = TET1;
        PSI2(i1:i1+2-1, i2:i2+8-1) = TET3;
        PSI1b(i1:i1+2-1, i2:i2+8-1) = TET1b;
        PSI2b(i1:i1+2-1, i2:i2+8-1) = TET3b;       
        PSI1t(i1:i1+2-1, i2t:i2t+6-1) = TET1t;
        PSI2t(i1t, i2t:i2t+6-1) = TET3t;         
    end
end
PSI1(end-2+1:end, 1:8) = TET1;
PSI2(end-2+1:end, 1:8) = TET3;
PSI1b(end-2+1:end, 1:8) = TET1b;
PSI2b(end-2+1:end, 1:8) = TET3b;
PSI1t(end-2+1:end, 1:6) = TET1t;
PSI2t(end, 1:6) = TET3t;

%
CS = zeros((2*q+1)*p,(8*q-2)*p);
CS(1:4*p,1:8*p) = [PSI1b;PSI2b];
CS(end-3*p+1:end,  end-6*p+1:end ) = [PSI1t;PSI2t];

for i = 2:q-1
    i1 = 2*p*(i-1)+1;
    i11 = 8*p*(i-1)+1;
    
    CS(i1:i1+2*p-1, i11:i11+8*p-1) = PSI1;
    
    i1 = 2*p*i+1;
    CS(i1:i1+2*p-1, i11:i11+8*p-1) = PSI2;
    
end


%CS(1:2*p, end-8*p+1:end) = PSI2;
CS = CS';

S = N*CS';

CSm1 = zeros(size(CS,1),size(CS,2));% CS-1
CSp1 = zeros(size(CS,1),size(CS,2));% CS+1
for i = 1:size(CS,1)
    for j = 1:size(CS,2)
        if CS(i,j) == -1
            CSm1(i,j) = -CS(i,j);
        elseif CS(i,j) == 1
            CSp1(i,j) = CS(i,j);
        end
    end
end

Sm1 = N*CSm1';
Sp1 = N*CSp1';

%PLOTTING STRINGS
for i = 1:size(S,2)
    plot3([Sm1(1,i),Sp1(1,i)], ...
        [Sm1(2,i),Sp1(2,i)], ...
        [Sm1(3,i),Sp1(3,i)],...
        'r','LineWidth',1.05)
end


xlabel('x', 'FontSize', 19);
ylabel('y', 'FontSize', 19);
zlabel('z', 'FontSize', 19);
set(gca,'FontSize',15);
%set(gca,'DefaultTextFontSize',24)
axis equal

n_s = size(CS,1)
%% Add the center point
% Find the points which z position = 0
sp_point = [];
for i = 1:size(N,2)
    if N(3,i) == 0
       sp_point = [sp_point , i] ;
    end
end


% Update N and Cb
N = [N center_node];

Cb = [CB zeros((2*q+1)*p,1)];
% Update Cs
CCC = - eye( 2*p*q+p );
size(CCC)
CCC = [ CCC  ones( (2*q+1)*p , 1)];
size(CCC)
% for i = 1:size(CCC,1)
%     if ismember(i,sp_point)
%         CCC(i,i) = 0;
%         CCC(i,size(CCC,2)) = 0;
%     end
% end
CCC(all(CCC == 0, 2),:) = [];
size(CS)
Cs = [CS  zeros(  size(CS,1) , 1 )];
size(Cs)
Cs = [Cs; CCC];


axis off
view(12,15)
%% Finding the number of nodes, strings and bars for N, Cs and Cb
nnodes = size(N,2);
n_s = size(Cs,1)
n_b = size(Cb,1);


%outputArg1 = inputArg1;
%outputArg2 = inputArg2;
end

