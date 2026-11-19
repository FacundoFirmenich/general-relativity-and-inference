function out = casimir_
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_A,par_B,par_C,par_D,par_W,par_d,par_n)
dydt=[kmrgd(2);
-par_A*kmrgd(2)-par_B*kmrgd(1)+par_D*kmrgd(3)-par_C*(kmrgd(1)+par_d)^par_n;
-par_W*kmrgd(4)+kmrgd(3)*(1-kmrgd(3)^2-kmrgd(4)^2);
par_W*kmrgd(3)+kmrgd(4)*(1-kmrgd(3)^2-kmrgd(4)^2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(casimir_);
y0=[0,0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_A,par_B,par_C,par_D,par_W,par_d,par_n)
jac=[ 0 , 1 , 0 , 0 ; - par_B - par_C*par_n*(kmrgd(1) + par_d)^(par_n - 1) , -par_A , par_D , 0 ; 0 , 0 , 1 - kmrgd(4)^2 - 3*kmrgd(3)^2 , - par_W - 2*kmrgd(3)*kmrgd(4) ; 0 , 0 , par_W - 2*kmrgd(3)*kmrgd(4) , 1 - 3*kmrgd(4)^2 - kmrgd(3)^2 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_A,par_B,par_C,par_D,par_W,par_d,par_n)
jacp=[ 0 , 0 , 0 , 0 , 0 , 0 , 0 ; -kmrgd(2) , -kmrgd(1) , -(kmrgd(1) + par_d)^par_n , kmrgd(3) , 0 , -par_C*par_n*(kmrgd(1) + par_d)^(par_n - 1) , -par_C*log(kmrgd(1) + par_d)*(kmrgd(1) + par_d)^par_n ; 0 , 0 , 0 , 0 , -kmrgd(4) , 0 , 0 ; 0 , 0 , 0 , 0 , kmrgd(3) , 0 , 0 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_A,par_B,par_C,par_D,par_W,par_d,par_n)
hess1=[ 0 , 0 , 0 , 0 ; -par_C*par_n*(kmrgd(1) + par_d)^(par_n - 2)*(par_n - 1) , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hess2=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hess3=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , -6*kmrgd(3) , -2*kmrgd(4) ; 0 , 0 , -2*kmrgd(4) , -2*kmrgd(3) ];
hess4=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , -2*kmrgd(4) , -2*kmrgd(3) ; 0 , 0 , -2*kmrgd(3) , -6*kmrgd(4) ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
hess(:,:,4) =hess4;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_A,par_B,par_C,par_D,par_W,par_d,par_n)
hessp1=[ 0 , 0 , 0 , 0 ; 0 , -1 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp2=[ 0 , 0 , 0 , 0 ; -1 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp3=[ 0 , 0 , 0 , 0 ; -par_n*(kmrgd(1) + par_d)^(par_n - 1) , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp4=[ 0 , 0 , 0 , 0 ; 0 , 0 , 1 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp5=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , -1 ; 0 , 0 , 1 , 0 ];
hessp6=[ 0 , 0 , 0 , 0 ; -par_C*par_n*(kmrgd(1) + par_d)^(par_n - 2)*(par_n - 1) , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp7=[ 0 , 0 , 0 , 0 ; - par_C*(kmrgd(1) + par_d)^(par_n - 1) - par_C*par_n*log(kmrgd(1) + par_d)*(kmrgd(1) + par_d)^(par_n - 1) , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_A,par_B,par_C,par_D,par_W,par_d,par_n)
tens31=[ 0 , 0 , 0 , 0 ; -par_C*par_n*(kmrgd(1) + par_d)^(par_n - 3)*(par_n - 1)*(par_n - 2) , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens32=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens33=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens34=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens35=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens36=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens37=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens39=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens310=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens311=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , -6 , 0 ; 0 , 0 , 0 , -2 ];
tens312=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , -2 ; 0 , 0 , -2 , 0 ];
tens313=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens314=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens315=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , -2 ; 0 , 0 , -2 , 0 ];
tens316=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , -2 , 0 ; 0 , 0 , 0 , -6 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,1,4) =tens34;
tens3(:,:,2,1) =tens35;
tens3(:,:,2,2) =tens36;
tens3(:,:,2,3) =tens37;
tens3(:,:,2,4) =tens38;
tens3(:,:,3,1) =tens39;
tens3(:,:,3,2) =tens310;
tens3(:,:,3,3) =tens311;
tens3(:,:,3,4) =tens312;
tens3(:,:,4,1) =tens313;
tens3(:,:,4,2) =tens314;
tens3(:,:,4,3) =tens315;
tens3(:,:,4,4) =tens316;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_A,par_B,par_C,par_D,par_W,par_d,par_n)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_A,par_B,par_C,par_D,par_W,par_d,par_n)
