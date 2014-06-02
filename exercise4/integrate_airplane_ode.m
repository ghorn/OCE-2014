function [xnext, Jx, Ju] = integrate_airplane_ode(x0,u)
    ts = 0.02;
    [dfdt,Jx_ct,Ju_ct] = airplane_ode(x0,u);
    xnext = x0 + dfdt*ts;
    Jx = eye(4) + ts*Jx_ct;
    Ju = ts*Ju_ct;
end


function [dfdt, Jx, Ja] = airplane_ode(xs,alpha)

params.ar = 10;
params.sref = 0.5;
params.mass = 2;

% px = xs(1);
% pz = xs(2);
vx = xs(3);
vz = xs(4);

% model params
ar = params.ar;
sref = params.sref;
mass = params.mass;
rho = 1.2;

% p = [px,pz];
v = [vx;vz];
ve = norm(v);

drag_e = -[vx;vz];
drag_e = drag_e/norm(drag_e);
lift_e = [vz; -vx];
lift_e = lift_e/norm(lift_e);

cla = 2*pi*10/12;
CL = cla*alpha;
CD = CL*CL/(pi*ar) + 0.01;

q = 1/2*rho*ve^2;
f = q*sref*(CL*lift_e + CD*drag_e);
fx = f(1);
fz = f(2) + 9.8*mass;

ddt_px = vx;
ddt_pz = vz;
ddt_vx = fx/mass;
ddt_vz = fz/mass;

dfdt = [ddt_px;
        ddt_pz;
        ddt_vx;
        ddt_vz];

psm = 1/2*rho*sref/mass;
Jx = [0,0,1,0;
      0,0,0,1;
      0,0,psm*(CL*vx*vz-CD*(2*vx*vx+vz*vz))/ve, psm*(CL*(vx*vx+2*vz*vz)-CD*vx*vz)/ve;
      0,0,-psm*(CD*vx*vz+CL*(2*vx*vx+vz*vz))/ve, -psm*(CL*vx*vz+CD*(vx*vx+2*vz*vz))/ve];

Ja = [0;
      0;
      psm*cla*ve*(vz-2*cla*vx*alpha/(pi*ar));
      -psm*cla*ve*(pi*ar*vx+2*cla*vz*alpha)/(pi*ar)];
     
return

end
