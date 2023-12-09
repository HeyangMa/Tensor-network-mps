%穷举  H=\sum -J{ij}*Si*Sj
close all
clear
clc
%% setting
nx=4;
ny=4;
J2=[];
% 1： 铁磁
J2(:,1)=textread('hxx.txt');
J2(:,2)=textread('hyy.txt');
h = -0.221;
beta_list=0.53;
pbc=0;
fbc=1;
%% main
zsum=zeros(1,length(beta_list));
for t=1:length(beta_list)
    beta=beta_list(t);
    %数组大小
    ising=-ones(1,nx*ny);
    nbor=neighbor_pbc(nx,ny);
    energy=zeros(2^(nx*ny),1);
    m=zeros(2^(nx*ny),1);
    %初始全为-1,计算
    [energy(1)]=calculate(ising,nbor,nx,ny,J2,h,pbc,fbc);
    z(1)=exp(-beta*energy(1));
    %2^n*n-1个值
    for i=1:length(energy)-1
        a=str2double(dec2bin(i));
        a=num2str(a,append('%0',num2str(nx*ny),'d'));
        for j=1:nx*ny
            if str2double(a(j)) == 1
                ising(j)=1;
            else
                ising(j)=-1;
            end
        end
        [energy(i+1)]=calculate(ising,nbor,nx,ny,J2,h,pbc,fbc);
        z(i+1)=exp(-beta*energy(i+1));
    end
    zsum(t)=sum(z);
    fprintf('pbc = %f      ',pbc);
    fprintf('fbc = %f      ',fbc);
    fprintf('Lx = %f      ',nx);
    fprintf('Ly = %f      ',ny);
    fprintf('beta = %f      ',beta);
    %fprintf('J = %f      ',J2);
    %fprintf('h = %f      ',h);
    fprintf('z = %f      ',zsum(t));
    fprintf('log(z) = %f\n',log(zsum(t)));
end
% figure(1);hold on;plot(beta_list,zsum,'k*-');set(gca,'Yscale','log');
%% generate neighbor
function [nbor]=neighbor_pbc(nx,ny)
nbor=zeros(nx*ny,4);
for ispin=1:nx*ny
    iy=fix((ispin-1)/nx)+1;
    ix=ispin-(iy-1)*nx;
    ixp=ix+1-fix(ix/nx)*nx;       %右侧点, x坐标
    iyp=iy+1-fix(iy/ny)*ny;       %上侧点, y坐标
    ixm=ix-1+fix((nx-ix+1)/nx)*nx; %左侧点, x坐标
    iym=iy-1+fix((ny-iy+1)/ny)*ny; %下侧点, y坐标
    nbor(ispin,1)=(iy-1)*nx+ixp;%右邻居
    nbor(ispin,2)=(iyp-1)*nx+ix;%上邻居
    nbor(ispin,3)=(iy-1)*nx+ixm;%左邻居
    nbor(ispin,4)=(iym-1)*nx+ix;%下邻居
end
end
%% calculate_energy
function [energy]=calculate(ising,nbor,nx,ny,J2,h,pbc,fbc)
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
    for i=1:nx*ny
        iy=fix((i-1)/nx)+1;
        ix=i-(iy-1)*nx;
        if ix~=nx && iy~=ny
            energy=energy - J2(i,1) * ising(i) * ising(nbor(i,1)) - J2(i,2) * ising(i) * ising(nbor(i,2));
        end
        if ix==nx && iy~=ny
            energy=energy - J2(i,2) * ising(i) * ( ising(nbor(i,2)) );
        end
        if iy==ny && ix~=nx
            energy=energy - J2(i,1) * ising(i) * ( ising(nbor(i,1)) );
        end
    end
    for i=1:nx*ny
        energy=energy - h * ising(i);
    end
end
end