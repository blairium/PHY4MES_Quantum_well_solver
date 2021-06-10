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
epsr=12.5; %Relative permittivity
m=0.067*9.1e-31; %Effective mass [kg]
k=8.617e-5; %Boltzmann constant [eV/K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculation parameters-Change these as necessary
a=3e-11; %mesh size [m] (Default is 30 pm)
mu=1   ; %Chemical potential [eV]
T=293    ; %Temperature [K]
Nout=100   ; %Outer passivating layer
Nin=66   ; %Inner passivating layer-seperates the two wells
Nw1=100   ; %First (left) well
Nw2=100   ; %Second (right) well
E_m=0.33 ;
E_out=0.33;
Vg=0; %Gate potential (not used for this calculation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%More parameters needed for calculations. Do not change
t0=(hbar^2)/(2*m*(a^2)*q); %Scaling factor
e0=q*a/eps0; %Scaling factor
kT=k*T; 
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
end;
ns=1e-4*sum(sum(n.*[zeros(Nout,1);ones(Nw1,1);zeros(Nin,1);ones(Nw2,1);zeros(Nout,1)]));
nn=1e-6*n./a;

for x=1:1:20
        Occ_x=log(1+exp((mu-D(x))./kT));
        Ed_x=P(:,x)*Occ_x.*P(:,x);
        N_x=2*n0*Ed_x;
        density(x)=1e-4*sum(sum(N_x.*[ones(Np,1)])); 
    end
    
    figure(2)
    plot(density);

figure(3)
plot(XX,P(:,1),XX,P(:,2),XX,P(:,3))
