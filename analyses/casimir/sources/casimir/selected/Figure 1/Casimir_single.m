function output = Casimir_single(wm,m,g,n,d,Pc,area,xinit,Fext)
global cds
%Initialize global variable cds, which is necessary to extend
%simulations. Matlab complains, but there is no other way.

%Matlab is a bitch and I need to manually convert everything to a double
%for some fucking reason
wm = double(wm);
m = double(m);
g = double(g);
n = double(n);
d = double(d);
Pc = double(Pc);
area = double(area);
xinit = double(xinit);
Fex = double(Fext);

%wm = 2*pi*16.30e6; %Hz - mechanical resonator frequency
t0 = (2*pi)/wm; %dummy variable: time
x0_scal = 1e-9; %dummy variable: position
%m = 3.96e-14; %pg - effective mass 
%g = 2*pi*170; %Hz - gamma/damping coefficient 
%n = -3.2; %used for matcont purposes 
n_r = double(-n); %calculated for our setup by Matthijs
%d = 18; %distance between drums, scaled
d_r = double(d*x0_scal); %distance between drums
%Pc = 2175.904;%158.55209773; %Pa - Casimir pressure from BCS model
%area = pi*(0.5*11.3e-6)^2; %area of the drums
Fc = Pc*area*d_r^n_r; %Casimir force at 40 nm
w = wm; %Hz - drive frequency
W = 1.02*w*t0;
period = (2*pi)/(w*t0); %period

init; %initialize matcont
fun=casimir_();

D = ((t0^2)*Fex)/(x0_scal*m);
C = (Fc*(t0^2))/(m*x0_scal^(n_r+1));
A = g*t0;
B = (t0^2)*(wm^2);

%initial simulation for transients
  
%parameter values and initial condition
p=[A;B;C;D;W;d;n];
X0=[xinit;0;1;0]; %-0.4101
P0=num2cell(p);

%the final point.
[~,y1]=ode45(fun{2},[0 10000],X0,[],P0{:}); %Use this if you know the final point already

while y1(end,2) > 0
    [~,y1]=ode45(fun{2},[0 period*0.15],y1(end,:),[],P0{:});
end

opt=odeset('MaxStep',.01); %better to have smaller steps for the continuation
x1=y1(end,:); %take the final point
[t,y]=ode45(fun{2},[0 period*1.025],x1,opt,P0{:}); 

%continuation of the periodic orbit 
%we provide the simulated response (y) to the initializer
tolerance=1e-7;
ap=5; %the active parameter is chosen as W
[x0,v0]=initOrbLC(@casimir_,t,y,p,ap,20,4,tolerance);
%v0 = x0;

%First, try a very brief continuation to see the direction
opt=contset;
opt.MaxNumPoints = 5;%increase if the plot doesn't cover the frequency range yet
opt=contset(opt,'InitStepsize',1e-6);
opt=contset(opt,'MaxStepSize',0.005);
opt.Backward = 1; %Start from initial frequency > 2pi (above peak)
opt.Singularities =1; %singularity detection on
opt.IgnoreSingularity = [1 2 4]; %ignore everything except LPC
opt.TSearchOrder = -1;

[xlc1,~,~,~,~]=cont(@limitcycle,x0,v0,opt);

freqs = xlc1(end,:);
if freqs(end) > freqs(1)
    opt.Backward = 0;
end

opt.MaxNumPoints = 1000;
[xlc1,vlc1,slc1,hlc1,flc1]=cont(@limitcycle,x0,v0,opt);

freqs2 = xlc1(end,:);
while freqs2(end) > 1.0 && freqs2(end) < 100
    disp(freqs2(end) +" > 1.0, adding 1000 extra points.")
    %Do a continuation from the previous limit cycle
    [xlc1,vlc1,slc1,hlc1,flc1]=cont(xlc1,vlc1,slc1,hlc1,flc1,cds);
    
    freqs2 = xlc1(end,:);
    
end

amps2 = (max(xlc1(1:4:end-2,:)) - min(xlc1(1:4:end-2,:)))/2;

output = [amps2,freqs2];