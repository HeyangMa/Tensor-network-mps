%穷举  H=\sum -J2*SiSj
close all
clear
clc
%% 设置参数
n=10;
J2=1;
h=1;
beta_list=1;
pbc=0;
fbc=1;
%% main
zsum=zeros(1,length(beta_list));
for t=1:length(beta_list)
    beta=beta_list(t);
    %数组大小
    ising=-ones(1,n);
    nbor=neighbor_pbc(n);
    energy=zeros(2^(n),1);
    m=zeros(2^(n),1);
    %初始全为-1,计算
    [energy(1)]=calculate(ising,nbor,n,J2,h,pbc,fbc);
    z(1)=exp(-beta*energy(1));
    %2^n*n-1个值
    for i=1:length(energy)-1
        a=str2double(dec2bin(i));
        a=num2str(a,append('%0',num2str(n),'d'));
        for j=1:n
            if str2double(a(j)) == 1
                ising(j)=1;
            else
                ising(j)=-1;
            end
        end
        [energy(i+1)]=calculate(ising,nbor,n,J2,h,pbc,fbc);
        z(i+1)=exp(-beta*energy(i+1));
    end
    zsum(t)=sum(z);
    fprintf('L = %f      ',n);
    fprintf('beta = %f      ',beta);
    fprintf('z = %f\n',zsum(t));
end
% figure(1);hold on;plot(beta_list,zsum,'k*-');set(gca,'Yscale','log');
%% generate neighbor
function [nbor]=neighbor_pbc(n)
nbor=zeros(n,2);
for ispin=1:n
    nbor(ispin,1)=ispin+1;%右邻居
    nbor(ispin,2)=ispin-1;%左邻居
    if ispin==1
        nbor(ispin,2)=n;  %左邻居
    end
    if ispin==n
        nbor(ispin,1)=1;  %右邻居
    end
end
end
%% calculate_energy
function [energy]=calculate(ising,nbor,n,J2,h,pbc,fbc)
%pbc
if pbc==1
    energy=0.0;
    for i=1:n
        energy=energy - J2 * ising(i) * ( ising(nbor(i,1))  )- h * ising(i);
    end
end

%fbc
if fbc==1
    energy=0.0;
    for i=1:n-1
        energy=energy - J2 * ising(i) * ( ising(nbor(i,1))  );
    end
    for i=1:n
        energy=energy - h * ising(i) ;
    end
end
end