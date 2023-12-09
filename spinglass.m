%穷举  H=\sum -J2*SiSj
close all
clear
clc
%% 改变n,变为不同尺寸
n=2;
J2=[];
% 1： 铁磁
J2(:,1)=textread('hxx.txt');
J2(:,2)=textread('hyy.txt');
h = 1;
beta_list=1;
pbc=0;
fbc=1;
%% main
zsum=zeros(1,length(beta_list));
for t=1:length(beta_list)
    beta=beta_list(t);
    %数组大小
    ising=-ones(1,n*n);
    nbor=neighbor_pbc(n);
    energy=zeros(2^(n*n),1);
    m=zeros(2^(n*n),1);
    %初始全为-1,计算
    [energy(1)]=calculate(ising,nbor,n,J2,h,pbc,fbc);
    z(1)=exp(-beta*energy(1));
    %2^n*n-1个值
    for i=1:length(energy)-1
        a=str2double(dec2bin(i));
        a=num2str(a,append('%0',num2str(n^2),'d'));
        for j=1:n*n
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
    fprintf('pbc = %f      ',pbc);
    fprintf('fbc = %f      ',fbc);
    fprintf('L = %f      ',n);
    fprintf('beta = %f      ',beta);
    %fprintf('J = %f      ',J2);
    %fprintf('h = %f      ',h);
    fprintf('z = %f      ',zsum(t));
    fprintf('log(z) = %f\n',log(zsum(t)));
end
% figure(1);hold on;plot(beta_list,zsum,'k*-');set(gca,'Yscale','log');
%% generate neighbor
function [nbor]=neighbor_pbc(n)
nbor=zeros(n*n,4);
for ispin=1:n*n
    iy=fix((ispin-1)/n)+1;
    ix=ispin-(iy-1)*n;
    ixp=ix+1-fix(ix/n)*n;       %右侧点, x坐标
    iyp=iy+1-fix(iy/n)*n;       %上侧点, y坐标
    ixm=ix-1+fix((n-ix+1)/n)*n; %左侧点, x坐标
    iym=iy-1+fix((n-iy+1)/n)*n; %下侧点, y坐标
    nbor(ispin,1)=(iy-1)*n+ixp;%右邻居
    nbor(ispin,2)=(iyp-1)*n+ix;%上邻居
    nbor(ispin,3)=(iy-1)*n+ixm;%左邻居
    nbor(ispin,4)=(iym-1)*n+ix;%下邻居
end
end
%% calculate_energy
function [energy]=calculate(ising,nbor,n,J2,h,pbc,fbc)
%pbc
if pbc==1
    energy=0.0;
    for i=1:n*n
        energy = energy - J2 * ising(i) * ( ising(nbor(i,1)) + ising(nbor(i,2)) ) - h * ising(i);
    end
end

%fbc
if fbc==1
    energy=0.0;
    for i=1:n*n
        iy=fix((i-1)/n)+1;
        ix=i-(iy-1)*n;
        if ix~=n && iy~=n
            energy=energy - J2(i,1) * ising(i) * ising(nbor(i,1)) - J2(i,2) * ising(i) * ising(nbor(i,2));
        end
        if ix==n && iy~=n
            energy=energy - J2(i,2) * ising(i) * ( ising(nbor(i,2)) );
        end
        if iy==n && ix~=n
            energy=energy - J2(i,1) * ising(i) * ( ising(nbor(i,1)) );
        end
    end
    for i=1:n*n
        energy=energy - h * ising(i);
    end
end
end