%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script is modified from one created by Alex Schenk                  %
% It taakes the inputs below a returns a series of plots                  %
%The script also creates folders for the plots named after the respective
%variables
%
clc;clear all; close all;

a = 3e-11 ;%
mu = 1 ;% Fermi level in eV
Temp = 298; % Kelvin
Nout = 100; %mesh points
Nw = 100; %mesh points
E_out = 4 ;%eV;
num_eig = 5 ;% Number of eigen states to display
saveloc = 'C:\Users\blair\Documents\University\2021\PHY4MES\Computational Lab 2 (2021)-20210504\figures\Q2'; %Save location
ED_eig_plots = 3; % Number of subplots of electron density, for some the density will be so low there's not much point plotting

%Loop through temp
for  Temp = [298, 373, 500]
    fprintf('Temperature of  %d K.\n',Temp)
    QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
    close all %Closes plot windows 
end
%Loop through well width
for Nout = [10, 50, 100, 200, 500] % Well dith of 1000 won't converge
    fprintf('Barrier Thickness of  %d nm.\n',Nout*.03)
    QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
    close all
end
% %Loop through barrier height
% for E_out = [1,2,3,4,5,8,10]
%     fprintf('Barrier height of  %d eV.\n',E_out)
%     QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
%     close all
% end
% %loop through fermi energies
% for mu = [0.001, 0.01, 0.1, 1]
%     fprintf('Fermi Energy of  %d eV.\n',mu)
%     QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
%     close all
% end
% %loop through barrier thickness
% for Nw = [10, 50, 100, 200, 500]
%     fprintf('Barrier Thickness of  %d nm.\n',Nw*.03)
%     disp('loop through barrier thickness')
%     QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
%     close all
% end
    
    
    
    

%[W,prob, E,nx_sum, U1, Ec, XX] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
function [W,prob, E,nx_sum, U1, Ec, XX] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots)


%%%%%%%%%%%%%%Inputs%%%%%%%%%%%%%%%%%%%%%
% a=3e-11; %mesh size [m] (Default is 30 pm) 3e-11;
% mu=  ; %Fermi level [eV]
% T=    ; %Temperature [K]
% Nout=   ; %Outer passivating layer
% Nw=   ; %Well width well
% E_out= ;
% Poisson-Schrodinger Solver    
% Iterative solver for semiconductor quantum well structures
% A Schenk 2015
% La Trobe University
% Altered from pre-existing code for PHY5PQA Matlab assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W - Wavefunction
%Prob - Probability densnity
% E - Energy Eigenvalue
% nx_sum - Electron density
% U1 - Potential
% Ec - 
% XX - Position

%Preset parameters-Do Not Change!
hbar=1.06e-34; %Plancks constant [Js]
q=1.6e-19; %Elementary charge [C]
eps0=8.85e-12; %Permitivity of free space [F/m]
epsr=4; %Relative permittivity
m=0.25*9.1e-31; %Effective mass [kg]
k=8.617e-5; %Boltzmann constant [eV/K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculation parameters-Change these as necessary
Vg=0; %Gate potential (not used for this calculation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%More parameters needed for calculations. Do not change
t0=(hbar^2)/(2*m*(a^2)*q); %Scaling factor
e0=q*a/eps0; %Scaling factor
kT=k*Temp; 
n0=m*kT*q/(2*pi*(hbar^2)); %2D DOS
Np=2*Nout+Nw; %layer thickness in units of mesh size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


XX=a*1e9*(1:1:Np);
Ec=[E_out*ones(Nout,1);0*ones(Nw,1);E_out*ones(Nout,1)];
T=(2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));
D2=epsr*((2*diag(ones(1,Np)))-(diag(ones(1,Np-1),1))-diag(ones(1,Np-1),-1));
iD2=inv(D2);
Vg=0;
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
ns=1e-4*sum(sum(n.*[zeros(Nout,1);ones(Nw,1);zeros(Nout,1)]));
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
 
    folder = strcat(saveloc,'\','Barrier_Height_',num2str(E_out),'\_eV_Well_width_', num2str(Nw) ,'\nm_Temp_' , num2str(Temp),'K\fermi_level', num2str(mu),'\barrier_thicknesss_',num2str(Nout),'\');
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
    

