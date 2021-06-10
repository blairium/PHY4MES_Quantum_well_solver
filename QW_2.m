%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script is modified from one created by Alex Schenk                  %
% It taakes the inputs below a returns a series of plots                  %
%The script also creates folders for the plots named after the respective
%variables
%
clc;clear all; close all;


a = 3e-11 ;% mesh point spacing
mu = 1 ;% Fermi level in eV
Temp = 298; % Kelvin
Nout = 100; %mesh points
Nin = 100; %Mesh point in inner barrier
Nw1 = 100; %mesh points
Nw2 = 100; %mesh points
E_m = 2 ; %eV
E_out = 4 ;%eV;
num_eig = 5 ;% Number of eigen states to display
saveloc = 'C:\Users\blair\Documents\University\2021\PHY4MES\Computational Lab 2 (2021)-20210504\figures\Q3'; %Save location
ED_eig_plots = 3; % Number of subplots of electron density, for some the density will be so low there's not much point plotting

%Nin=0
QW2( a, mu, Temp,Nout, 0, Nw1, Nw2, E_m, E_out, num_eig, saveloc, ED_eig_plots);

%E_m=E_out, Nw1=100, Nw2=170, Nin range 2nm to 20nm
for Nin = [66, 111, 222, 333, 444, 555, 666]
    QW2( a, mu, Temp,Nout, Nin, 100, 170, E_out, E_out, num_eig, saveloc, ED_eig_plots)
    close all
end

for E_m = [1,2,3,4,5]
    Nw1= 100;
    Nw2 = 200;
    QW2( a, mu, Temp,Nout, Nin, Nw1, Nw2, E_m, E_out, num_eig, saveloc, ED_eig_plots);
    close all
end

function [W,prob, E,nx_sum, U1, Ec, XX] = QW2(a, mu, Temp,Nout, Nin, Nw1, Nw2, E_m, E_out, num_eig, saveloc, ED_eig_plots)
%Input Calculation parameters-Change these as necessary
% a=3e-11; %mesh size [m] (Default is 30 pm)
% mu=   ; %Chemical potential [eV]
% T=    ; %Temperature [K]
% Nout=   ; %Outer passivating layer
% Nin=   ; %Inner passivating layer-seperates the two wells
% Nw1=   ; %First (left) well
% Nw2=   ; %Second (right) well
% E_m= ;Energy height of the internal barrier layer
% E_out= ;Energy height of the outer barrier
Vg=0; %Gate potential (not used for this calculation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W - Wavefunction
%Prob - Probability densnity
% E - Energy Eigenvalue
% nx_sum - Electron density
% U1 - Potential
% Ec - 
% XX - Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Poisson-Schrodinger Solver    
% Iterative solver for semiconductor quantum well structures
% A Schenk 2015
% La Trobe University
% Altered from pre-existing code for PHY5PQA Matlab assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preset parameters-Do Not Change!
hbar=1.06e-34; %Plancks constant [Js]
q=1.6e-19; %Elementary charge [C]
eps0=8.85e-12; %Permitivity of free space [F/m]
epsr=4; %Relative permittivity
m=0.25*9.1e-31; %Effective mass [kg]
k=8.617e-5; %Boltzmann constant [eV/K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%More parameters needed for calculations. Do not change
t0=(hbar^2)/(2*m*(a^2)*q); %Scaling factor
e0=q*a/eps0; %Scaling factor
kT=k*Temp; 
n0=m*kT*q/(2*pi*(hbar^2)); %2D DOS
Np=2*Nout+Nin+Nw1+Nw2; %layer thickness in units of mesh size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matrices needed for setting up calculation
XX=a*1e9*[1:1:Np]; %Position in nm
Ec=[E_out*ones(Nout,1);0*ones(Nw1,1);E_m*ones(Nin,1);0*ones(Nw2,1);E_out*ones(Nout,1)]; %Conduction band position
T=(2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1)); 
D2=epsr*((2*diag(ones(1,Np)))-(diag(ones(1,Np-1),1))-diag(ones(1,Np-1),-1));
iD2=inv(D2);

Ubdy=-4*[Vg;zeros(Np-2,1);Vg]; 
U0=iD2*Ubdy;
U1=1e-9*ones(Np,1);UU=U1;change=1;
while change>1e-6
    U1=U1+0.1*(UU-U1);
    [P,D]=eig(T+diag(Ec)+diag(U1));D=diag(D);
    rho=log(1+exp((mu-D)./kT)); rho=P*diag(rho)*P';
    n=2*n0*diag(rho);
    UU=U0+iD2*e0*n;
    change=max(max((abs(UU-U1))));
    U=Ec+U1;
    disp(change)
end
ns=1e-4*sum(sum(n.*[zeros(Nout,1);ones(Nw1,1);zeros(Nin,1);ones(Nw2,1);zeros(Nout,1)]));
nn=1e-6*n./a;
function w_n = w_n(x)
        w_n = P(:,x);
    end
    function prob = w2(x)
        prob = w_n(x).^2;
    end

    function EDx = ED(x) %Electron Density
        EDx = D(x);
    end
    function Occ_x = OC(x)
        Occ_x=log(1+exp((mu-D(x))./kT));
    end

    function Ed_x = EDx(x)
        Ed_x=P(:,x).^2 *OC(x) ;%*P(:,x)';
    end
        
    function N_x = NX(x)
        N_x=2*n0*EDx(x);
    end
        
    function nx_sum = NXSUM(x) %Sum electron density
        nx_sum=1e-4*sum(sum(NX(x).*[ones(Np,1)]));
    end


% shift xx so that interface falls at origin for all dopant conc
    XX = XX-min(XX);
    XX = XX-max(XX)/2;
 
    folder = strcat(saveloc,'\','Barrier_Height_',num2str(E_out),'\_eV_Well_1_', num2str(Nw1),'nm_','Well_2_',num2str(Nw2) ,'nm\Temp_' , num2str(Temp),'K\fermi_level', num2str(mu),'\barrier_thicknesss_',num2str(Nout),'\');
    mkdir(folder)
    figure(1)
    
    hold on
    for x = 1:num_eig
        label = strcat('Eigenstate number ', num2str(x)); %Labels each state
        plot(XX,w_n(x),'DisplayName',label,'LineWidth',2)
    end

    title('Barrier Height','interpreter','latex')
    xlabel(' Position($n m$)','interpreter','latex')
    ylabel('Energy eV','interpreter','latex')


    legend('Location','southeast')
    saveas(gcf,strcat(folder,'Wavefunctions.png'))

    hold off
    

    figure(2)
    hold on
    for x = 1:num_eig
        label = strcat('Eigenstate number ', num2str(x));
        plot(XX,w2(x),'DisplayName',label,'LineWidth',2)
    end

    title('Probability Density','interpreter','latex')
    xlabel(' Position($n m$)','interpreter','latex')
    ylabel('Energy eV','interpreter','latex')


    legend('Location','southeast')
    saveas(gcf,strcat(folder,'prob_density.png'))

    hold off
    

    figure(3)
    hold on
    for x = 1:ED_eig_plots
        subplot(ED_eig_plots,1,x)
        plot(XX,EDx(x),'DisplayName',label,'LineWidth',2)
        title(strcat('Electron Density for eigenvalue ',' ',num2str(x)),'interpreter','latex');
        xlabel(' Position($n m$)','interpreter','latex')
        ylabel('Energy eV','interpreter','latex')

    end



    %legend('Location','southeast')
    saveas(gcf,strcat(folder,'elec_density.png'))

    hold off

    figure(4)
    hold on

    plot(XX,U1,'LineWidth',2)

    title('Potential vs location','interpreter','latex')
    xlabel('Position($n m$)','interpreter','latex')
    ylabel('Energy eV','interpreter','latex')


    legend('Location','southeast')
    saveas(gcf,strcat(folder,'potential.png'))

    hold off  


    figure(5)
    hold on

    plot(XX,Ec,'DisplayName','Initial Potential Profile','LineWidth',2)
    plot(XX,U,'DisplayName','Final Potential Profile','LineWidth',2)

    title('Sum Electron Density ','interpreter','latex')
    xlabel(' Position($n m$)','interpreter','latex')
    ylabel('Energy eV','interpreter','latex')


    legend('Location','southeast')
    saveas(gcf,strcat(folder,'profile.png'))

    hold off

end

