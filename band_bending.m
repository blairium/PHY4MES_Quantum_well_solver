% Band Bending Calculation      
% Poisson/Drift solver for semiconductor heterojunctions 
% A Schenk 2015
% La Trobe University


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All charge carrier and dopant concentrations are defined per cm^3

%Important Constants
q     = 1.602E-19;        % Elementary Charge
kb    = 1.381E-23;         % Boltzmann constant [J/K]
eps   = 8.85e-14;         % Permittivity of free space [F/cm]
T     = 300;              % Temperature [K]
Vt    = kb*T/q;           % "Thermal Voltage" [eV]

%N type side
ni_N=   ;      %Intrinsic carrier concentration on N doped side, [per cm^3]
Nd=   ;                %Donor concentration, [per cm^3]
Nc=   ;            %Conduction band DOS, [per cm^3]
dEc=Vt*log(Nc/Nd);%Position of Fermi level relative to conduction band, [eV]
chi_N=   ;       %Electron Affinity of N doped side, used for band alignment, [eV]
Egap_N=   ;        % Band Gap of N doped Side [eV]
kappa_N=   ;      %Relative Permittivity of N doped side

%P type side
ni_P=   ;       %Intrinsic carrier concentration on P doped side [per cm^3]
Na=    ;                %Acceptor concentration [per cm^3]
Nv=     ;            %Valence Band DOS [per cm^3]
dEv=Vt*log(Nv/Na);  %Position of Fermi level relative to conduction band [eV]
chi_P=   ;          %Electron Affinity of P doped side, used for band alignment [eV]
Egap_P=   ;        %Band Gap of P doped side [eV]
kappa_P=   ;      %Relative Permittivity of P doped side   


Vbi=(chi_P-chi_N)+Egap_P-dEc-dEv;  %Calculates difference in position between the Fermi levels of P and N doped sides (Contact Potential Difference) [eV]
V_Jp=Vbi*kappa_N*Nd/(kappa_N*Nd+kappa_P*Na); %Weighted potential in P type material, used as starting parameter [eV]
V_Jn=V_Jp*kappa_P*Na/(kappa_N*Nd); %Same as above for N type [eV]
Wp=(2*kappa_P*eps*V_Jp/(q*Na))^0.5; %Depletion width in P type material based on Depletion approximation-used as starting parameter [cm]
Wn=(2*kappa_N*eps*V_Jp/(q*Nd))^0.5; %Depletion width in N type material based on Depletion approximation-used as starting parameter [cm]
Ldn = sqrt(eps*kappa_N*Vt/(q*Nd)); %Debye length of electric field in N side, used as a scaling factor during calculations [cm]
Ldp = sqrt(eps*kappa_P*Vt/(q*Na)); %Debye length of electric field in P side, used as a scaling factor during calculations [cm]
Ldin = sqrt(eps*kappa_N*Vt/(q*ni_N)); %Debye length of intrinsicly doped material on N side [cm]
Ldip= sqrt(eps*kappa_P*Vt/(q*ni_P)); %Debye length of intrinsicly doped material on P side [cm]

%Defines the intrinsic Debye length to use for scaling based on which is
%larger, in order to perform a sensible calculation
Ldi = 0;
if(Ldi < Ldin)
    Ldi = Ldin;
end
if(Ldi < Ldip)
    Ldi=Ldip;
end

%Defines length scale of system based on Depletion approximation parameters

x_max = 0;
if(x_max < Wn)
    x_max = Wn;
end
if(x_max < Wp)
    x_max = Wp;
end
x_max = 100*x_max;       %This value is simply a "best guess"; you may need to increase beyond 100*x_max if you cannot get a flat region in your band

% Setting the grid size based on the extrinsic Debye lengths (effectively decay length of the electric field) %

dx = Ldn;
if(dx > Ldp)
    dx=Ldp;
end
dx = dx/50;   % Step size of the simulation. Scaling factor of 1/50 can be adjusted to speed up the calculation (increasing step size), at the cost of resolution

% Calculate the required number of grid points and renormalize dx %

n_max = x_max/dx;
n_max = round(n_max);

dx = dx/Ldi;    % Renormalize lengths with Ldi


%Determines which intrinsic carrier concentration to use for scaling,
%based on which is bigger, in order to produce sensible numbers
nscale=ni_N
if (ni_P>ni_N)
    nscale=ni_P
end

% Set up the doping C(x)=Nd(x)-Na(x) that is normalized with ni %
% Initial condition for holes/electrons %

for i = 1:n_max
    if(i <= n_max/2)
        dop(i) = -Na/nscale;
    elseif(i > n_max/2)
        dop(i) = Nd/nscale;
    end
end

% Initialize the potential based on the requirement of charge
% neutrality throughout the whole structure
format('long')
for i = 1: n_max
    zz = 0.5*dop(i);
    fi(i) = log(abs(zz));
    n(i) = abs(zz);
    p(i) = 1/abs(zz);
end

delta_acc = 1E-5;               % Preset the Tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                      %%
%%               EQUILIBRIUM  SOLUTION PART BEGINS                      %%
%%                                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(A) Define the elements of the coefficient matrix for the internal nodes and
%    initialize the forcing function

dx2 = dx*dx;
for i = 1: n_max
    a(i) = 1/dx2;
    c(i) = 1/dx2;
    b(i) = -(2/dx2+exp(fi(i))+exp(-fi(i)));
    f(i) = exp(fi(i)) - exp(-fi(i)) - dop(i) - fi(i)*(exp(fi(i))+exp(-fi(i)));
end


%(B) Define the elements of the coefficient matrix and initialize the forcing
%    function at the ohmic contacts 

a(1) = 0;
c(1) = 0;
b(1) = 1;
f(1) = fi(1);
a(n_max) = 0;
c(n_max) = 0;
b(n_max) = 1;
f(n_max) = fi(n_max);

%(C)  Start the iterative procedure for the solution of the linearized Poisson
%     equation using LU decomposition method:

flag_conv = 0;		           % convergence of the Poisson loop
k_iter= 0;
while(~flag_conv)            
    k_iter = k_iter + 1; 
    
    alpha(1) = b(1);
    for i=2:n_max
        beta(i)=a(i)/alpha(i-1);
        alpha(i)=b(i)-beta(i)*c(i-1);
    end
    
% Solution of Lv = f %    

    v(1) = f(1);
    for i = 2:n_max
        v(i) = f(i) - beta(i)*v(i-1);
    end
     
% Solution of U*fi = v %    

    temp = v(n_max)/alpha(n_max);
    delta(n_max) = temp - fi(n_max);
    fi(n_max)=temp;
    for i = (n_max-1):-1:1       %delta%
        temp = (v(i)-c(i)*fi(i+1))/alpha(i);
        delta(i) = temp - fi(i);
        fi(i) = temp;
    end
    
    delta_max = 0;
    
    for i = 1: n_max
        xx = abs(delta(i));
        if(xx > delta_max)
            delta_max=xx;
        end
    end

   
% Test convergence and recalculate forcing function and 
% central coefficient b if necessary
    
    if(delta_max < delta_acc)
        flag_conv = 1;
    else
        for i = 2: n_max-1
            b(i) = -(2/dx2 + exp(fi(i)) + exp(-fi(i)));
            f(i) = exp(fi(i)) - exp(-fi(i)) - dop(i) - fi(i)*(exp(fi(i)) + exp(-fi(i)));
        end
    end
end

xx1(1) = dx*1e4;

%This loop calculates most of the values of interest in this assignment

for i = 2:n_max-1 
    if(i <= n_max/2)
    Ec(i) = -chi_P+dEc+dEv- Vt*fi(i); %Conduction band position [eV]
    Ev(i)=Ec(i)-Egap_P; %Valence band position [eV]
    ro(i) = -ni_P*(exp(fi(i)) - exp(-fi(i)) - dop(i)); %Charge density [C/cm^3]
    else  Ec(i) = -chi_N+dEv+dEc - Vt*fi(i); %Conduction band position [eV]
    Ev(i)=Ec(i)-Egap_N; %Valence band position [eV]
        ro(i) = -ni_N*(exp(fi(i)) - exp(-fi(i)) - dop(i)); %Charge density [C/cm^3]
    end
    % This code calculates the electric field using two finite element
    % approximations. Ideally the results should be identical; if these two
    % quantities are significantly different then there may be an issue
    % with your calculation
    el_field1(i) = -(fi(i+1) - fi(i))*Vt/(dx*Ldi); %Electric field, approx 1
    el_field2(i) = -(fi(i+1) - fi(i-1))*Vt/(2*dx*Ldi); %Electric field, approx 2
    
    n(i) = exp(fi(i));  %scaled Electron density
    p(i) = exp(-fi(i)); %scaled Hole density
    xx1(i) = xx1(i-1) + dx*Ldi*1e4; %Calculating x position
end


%The following lines fix missing data points at the start/end of each
%dataset, which arise due to the finite element nature of the calculation
Ec(1) = Ec(2);
Ec(n_max) = Ec(n_max-1);
Ev(1)=Ev(2);
Ev(n_max)=Ev(n_max-1);
xx1(n_max) = xx1(n_max-1) + dx*Ldi*1e4;
el_field1(1) = el_field1(2);
el_field2(1) = el_field2(2);
el_field1(n_max) = el_field1(n_max-1);
el_field2(n_max) = el_field2(n_max-1);
ro(1) = ro(2);
ro(n_max) = ro(n_max-1);

% This loop a previous rescaling of the carrier density 
nf=n*nscale; %Final electron density
pf=p*nscale; %Final hole density

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%     Add your plotting codes after this line
% Variables which may be of interest for plotting are:
% Ec-Conduction band position [eV]
% Ev-Valence band position [eV]
% nf-electron density [per cm^3]
% pf-hole density [per cm^3]
% ro-total charge density [C/cm^3)
% el_field1-electric field, approximation 1 [V/cm]
% el_field2-electric field, approximation 2 [V/cm]
% xx1 should be used as your distance coordinates. The scaling in this code
% means that this is measured in micrometers
% 
