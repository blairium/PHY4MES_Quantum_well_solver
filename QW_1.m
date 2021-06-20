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
Nw = 412; %mesh points
E_out = 4 ;%eV;
num_eig = 1 ;% Number of eigen states to display
saveloc = 'C:\Users\blair\Documents\University\2021\PHY4MES\Computational Lab 2 (2021)-20210504\figures'; %Save location
ED_eig_plots = 3; % Number of subplots of electron density, for some the density will be so low there's not much point plotting
folder = strcat(saveloc,'\Q2\');
mkdir(folder)


hold on
for    Temp = 100:10:1000
    fprintf('Temp  %d K.\n',Temp)
    [~,~,NX_SUM1,NX_SUM2, NX_SUM3, NX_SUM4, ~, ~, ~, ~] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots, plot);
    scatter(Temp,NX_SUM1,'filled', 'blue');
    scatter(Temp,NX_SUM2,'filled', 'red');
    scatter(Temp,NX_SUM3,'filled', 'green');
    scatter(Temp,NX_SUM4,'filled', 'yellow');
    title('Sum Electron Density vs Barrier Thickness ','interpreter','latex')
    xlabel('Inner Barrier Thickness($nm$)','interpreter','latex')
    ylabel('Electron Density','interpreter','latex')
    legend('First Eigenstate','Second Eigenstate', 'Third Eigenstate', 'Fourth Eigenstate','location','best')

    
    
end
saveas(gcf,strcat(folder,'nx_sum_vs_inner_barrier.png'))
hold off
% figure(1)
% hold on
% for  Temp = 100:10:1000 %Temperatures less that 100ish cause problems, return inf during convergence loop
% 
%     [W,prob,nx_sum, U1, Ec, XX, EDxx] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
%     yyaxis left
%     scatter(Temp,nx_sum,'filled','blue');
%     title({'Sum Electron Density and Energy Eigenvalues';'vs Temperature '},'interpreter','latex')
%     xlabel(' Temperature ($^\circ$K)','interpreter','latex')
%     ylabel('Electron Density','interpreter','latex')
%     
%     yyaxis right
%     scatter(Temp,EDxx,'filled', 'red');
%     ylabel('Electron Energy Eigenvalue eV','interpreter','latex')
%     saveas(gcf,strcat(folder,'nx_sum_vs_temp.png'))
%     
% end
% hold off

% figure(2)
% hold on
% for  mu = 0:0.01:2
% 
%     [W,prob,nx_sum, U1, Ec, XX, EDxx] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
%     yyaxis left
%     scatter(mu,nx_sum,'filled','blue');
%     title({'Sum Electron Density and Energy Eigenvalues';'vs Fermi level '},'interpreter','latex')
%     xlabel(' Fermi Level ($eV$)','interpreter','latex')
%     ylabel('Electron Density','interpreter','latex')
%     
%     yyaxis right
%     scatter(mu,EDxx,'filled', 'red');
%     ylabel('Electron Energy Eigenvalue eV','interpreter','latex')
%     saveas(gcf,strcat(folder,'nx_sum_vs_fermi.png'))
%     
% end
% hold off
% 
% 
% figure(3)
% hold on
% for  Nw = 10:1:300
% 
%     [W,prob,nx_sum, U1, Ec, XX, EDxx] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
%     yyaxis left
%     scatter(Nw*0.03,nx_sum,'filled','blue');
%     title({'Sum Electron Density and Energy Eigenvalues';'vs well width'},'interpreter','latex')
%     xlabel('Well width ($nm$)','interpreter','latex')
%     ylabel('Electron Density','interpreter','latex')
%     
%     yyaxis right
%     scatter(Nw*0.03,EDxx,'filled', 'red');
%     ylabel('Electron Energy Eigenvalue eV','interpreter','latex')
%     saveas(gcf,strcat(folder,'nx_sum_vs_well.png'))
%     
% end
% hold off

%figure(4)
% hold on
% for  Nout = 50:5:300
% 
%     [W,prob,nx_sum, U1, Ec, XX, EDxx] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
%     yyaxis left
%     scatter(Nw*0.03,nx_sum,'filled','blue');
%     title({'Sum Electron Density and Energy Eigenvalues';'vs barrier thickness'},'interpreter','latex')
%     xlabel('Well width ($nm$)','interpreter','latex')
%     ylabel('Electron Density','interpreter','latex')
%     
%     yyaxis right
%     scatter(Nw*0.03,EDxx,'filled', 'red');
%     ylabel('Electron Energy Eigenvalue eV','interpreter','latex')
%     saveas(gcf,strcat(folder,'nx_sum_vs_thick_barrier.png'))
%     
% end
% hold off

% %Loop through well width
% for Nout = [10, 50, 100, 200, 500] % Well dith of 1000 won't converge
%     fprintf('Barrier Thickness of  %d nm.\n',Nout*.03)
%     QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots);
%     close all
% end
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





    



%[W,prob, E,nx_sum, U1, Ec, XX, EDX] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots, n);
function [W,prob,NX_SUM1,NX_SUM2, NX_SUM3, NX_SUM4, U1, Ec, XX, EDxx] = QW1(a, mu, Temp, Nout, Nw, E_out, num_eig, saveloc, ED_eig_plots, plot)


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
	W = w_n(num_eig);
    function prob = w2(x)
        prob = w_n(x).^2;
    end
    prob = w2(num_eig);
    function EDx = ED(x) %Electron Density
        EDx = D(x);
    end
    EDxx = ED(num_eig);
    function Occ_x = OC(x)
        Occ_x=log(1+exp((mu-D(x))./kT));
    end

    function Ed_x = EDx(x)
        Ed_x=P(:,x).^2 *OC(x) ;%*P(:,x)';
    end
    EDX = EDx(num_eig);
    function N_x = NX(x)
        N_x=2*n0*EDx(x);
    end
        
    function nx_sum = NXSUM(x) %Sum electron density
        nx_sum=1e-4*sum(sum(NX(x).*[ones(Np,1)]));
    end
        
    NX_SUM1 = NXSUM(1);
    NX_SUM2 = NXSUM(2);
    NX_SUM3 = NXSUM(3);
    NX_SUM4 = NXSUM(4);
   


% shift xx so that interface falls at origin for all dopant conc
    XX = XX-min(XX);
    XX = XX-max(XX)/2;
   
 
%     folder = strcat(saveloc,'\','Barrier_Height_',num2str(E_out),'\_eV_Well_width_', num2str(Nw) ,'\nm_Temp_' , num2str(Temp),'K\fermi_level', num2str(mu),'\barrier_thicknesss_',num2str(Nout),'\');
%     mkdir(folder)
%     figure(1)
% 
%     hold on
%     for x = 1:num_eig
%         label = strcat('Eigenstate number ', num2str(x)); %Labels each state
%         plot(XX,w_n(x),'DisplayName',label,'LineWidth',2)
%     end
% 
%     title('Barrier Height','interpreter','latex')
%     xlabel(' Position($n m$)','interpreter','latex')
%     ylabel('Energy eV','interpreter','latex')
% 
% 
%     legend('Location','southeast')
%     saveas(gcf,strcat(folder,'Wavefunctions.png'))
% 
%     hold off
% 
% 
%     figure(2)
%     hold on
%     for x = 1:num_eig
%         label = strcat('Eigenstate number ', num2str(x));
%         plot(XX,w2(x),'DisplayName',label,'LineWidth',2)
%     end
% 
%     title('Probability Density','interpreter','latex')
%     xlabel(' Position($n m$)','interpreter','latex')
%     ylabel('Energy eV','interpreter','latex')
% 
% 
%     legend('Location','southeast')
%     saveas(gcf,strcat(folder,'prob_density.png'))
% 
%     hold off
% 
% 
%     figure(3)
%     hold on
%     for x = 1:ED_eig_plots
%         subplot(ED_eig_plots,1,x)
%         plot(XX,EDx(x),'DisplayName',label,'LineWidth',2)
%         title(strcat('Electron Density for eigenvalue ',' ',num2str(x)),'interpreter','latex');
%         xlabel(' Position($n m$)','interpreter','latex')
%         ylabel('Energy eV','interpreter','latex')
% 
%     end
% 
% 
% 
%     %legend('Location','southeast')
%     saveas(gcf,strcat(folder,'elec_density.png'))
% 
%     hold off
% 
%     figure(4)
%     hold on
% 
%     plot(XX,U1,'LineWidth',2)
% 
%     title('Potential vs location','interpreter','latex')
%     xlabel('Position($n m$)','interpreter','latex')
%     ylabel('Energy eV','interpreter','latex')
% 
% 
%     legend('Location','southeast')
%     saveas(gcf,strcat(folder,'potential.png'))
% 
%     hold off  
% 
% 
%     figure(5)
%     hold on
% 
%     plot(XX,Ec,'DisplayName','Initial Potential Profile','LineWidth',2)
%     plot(XX,U,'DisplayName','Final Potential Profile','LineWidth',2)
% 
%     title('Sum Electron Density ','interpreter','latex')
%     xlabel(' Position($n m$)','interpreter','latex')
%     ylabel('Energy eV','interpreter','latex')
% 
% 
%     legend('Location','southeast')
%     saveas(gcf,strcat(folder,'profile.png'))
% 
%     hold off


end
    

