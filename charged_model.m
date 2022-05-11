%% Ellipsoid
a = 1;
b = 1;
c = 3;
% arbitrary set time
dt = 0.01;
t = 1:dt:25;
n = length(t);

% frequency variable longitudinally around ellipsoid
omega = 0.0;
% frequency variable latitudinally around ellipsoid
delta = 0.25;

% set equal spacing using 30 points
theta = linspace(0,pi,30);
phi = linspace(0,2*pi,30);

[th, ph] = meshgrid(theta, phi);

% allocate memory
% ellipsoid surface
x = zeros(length(th),length(ph));
y = zeros(length(th),length(ph));
z = zeros(length(th),length(ph));

% particle path
X = zeros(length(th),length(ph));
Y = zeros(length(th),length(ph));
Z = zeros(length(th),length(ph));

% Plot the surface (polar reference)
figure(1);
x = a.*sin(th).*cos(ph);
y = b.*sin(th).*sin(ph);
z = c.*cos(th);

ellip = surf(x,y,z);
hold on

% from center of ellipsoid
% This path shows a particle starting at the bottom pole of the ellipsoid 
% then travelling upward to the top pole of the ellipsoid 
theta = linspace(-pi/2,pi/2,n);

% theta = linspace(pi/2,-pi/2,n); % travelling from top pole to bottom pole
% theta = delta*t;

% path of particle (respect to center)
X = a.*cos(theta).*cos(omega*t);
Y = b.*cos(theta).*sin(omega*t);
Z = c.*sin(theta);

plot3(X,Y,Z,'r','linewidth',2);
% This poroportionally scales the 3 dimensional axes so the ellipsoid appears
% properly elongated
h = get(gca,'DataAspectRatio');
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
      set(gca,'DataAspectRatio',[1 1 h(3)])
end
hold off

% % velocity components
Vx = -omega.*a.*cos(theta).*sin(omega*t);
Vy = omega.*b.*cos(theta).*cos(omega*t);
Vz = c.*cos(theta);
% acceleration components
Ax = -(omega^2).*a.*cos(theta).*cos(omega*t);
Ay = -(omega^2).*b.*cos(theta).*sin(omega*t);
Az = c.*cos(theta);

%////////////////////////////////////////////
% These describe the particle elements when theta depends on time
% theta velocity components
% Vx = -omega.*a.*cos(delta*t).*sin(omega*t) - delta.*a.*sin(delta*t).*cos(omega*t);
% Vy = omega.*b.*cos(delta*t).*cos(omega*t) - delta.*b.*sin(delta*t).*sin(omega*t);
% Vz = -delta.*c.*sin(delta*t);
% % theta acceleration components
% Ax = (omega.*delta.*a.*sin(delta*t).*sin(omega*t) - (omega^2).*a.*cos(delta*t).*cos(omega*t)) + (-(delta^2).*a.*cos(delta*t).*cos(omega*t) + delta.*omega.*a.*sin(delta*t).*cos(omega*t));
% Ay = (-omega.*delta.*b.*sin(delta*t).*cos(omega*t) - (omega^2).*b.*cos(delta*t).*sin(omega*t)) + (-(delta^2).*b.*cos(delta*t).*sin(omega*t) + delta.*omega.*b.*cos(delta*t).*sin(omega*t));
% Az = -(delta^2).*c.*cos(delta*t);

% Observer at center of shape
figure(6);
hold on
r = zeros(3,n);
r(1,1) = 2;
scatter3(r(1,:),r(2,:),r(3,:),'o');
hold off

% position of particle
w = [X; Y; Z];
v = [Vx; Vy; Vz];
a = [Ax; Ay; Az];

% R is vector between observer(r) and particle position(w)
R = r - w;
% figure(4);
% plot(phi,R);

q = 1; % charge of particle
speed_light = 1E8; %3.0E8; % m/s
eps0 = 1; %8.8541878E-12; % C^2/Nm^2
k = q/(4*pi()*eps0);

u = (R./norm(R)).*speed_light - v;

E = zeros(3,n);

% calculation of electric field
for i=1:n
    E(:,i) = (k*norm(R(:,i))./(dot(R(:,i),u(:,i)).^3)) .* ((speed_light^2 - norm(v(:,i))^2).*u(:,i) + cross( R(:,i), (cross(u(:,i),a(:,i))) ));
end
% calculation of magnetic field
B = zeros(3,n);
B = 1/c^2 .* cross(v,E);

figure(2);
hold on
title('Electric Field');
xlabel('Time');
ylabel('Electric Field Magnitude');
% Z field
plot(t,E(1,:));
% Y field
plot(t,E(2,:));
% Z field
plot(t,E(3,:));
legend('E_x','E_y','E_z');
hold off

figure(3);
hold on
title('Magnetic Field');
xlabel('Time');
ylabel('Magnetic Field Magnitude')
% X field
plot(t,B(1,:));
% Y field
plot(t,B(2,:));
% Z field
plot(t,B(3,:));
legend('B_x','B_y','B_z');
hold off

%% Poynting Vector
permiability = 1E-6; %mkg/sA^2
%poynting vector calculation
s = zeros(3,n);
s = (1/permiability).*cross(E,B);
figure(4)
hold on
title('Poynting Vector');
xlabel('Time');
ylabel('Poynting Vector Magnitude');
plot(t,s(1,:));
plot(t,s(2,:));
plot(t,s(3,:));
legend('S_x','S_y','S_z');

%% Radiated Power eq. 11.7
P = zeros(1,length(t));
for i=1:length(t)
    P(i) = (permiability*q^2*norm(a(:,i))^2/(6*pi()*speed_light));
end
figure(5)
hold on
title('Power');
xlabel('Time');
ylabel('Power');
plot(t,P);

%% Electrostatics --> Electrodynamics
q1 = 1;
% new charge location of observer at t = 0
r_charge = zeros(3,n);
r_charge(1,1) = 2;

q2 = 0.001;
vel2 = zeros(3,n);
acc = zeros(3,n);
m2 = 0.001;

% Using Lorentz eq
for i=1:n
    % get new distance between observer and charged particle
    R(:,i) = r_charge(:,i) - w(:,i);
    % unit vector calculation
    u(:,i) = (R(:,i)./norm(R(:,i))).*speed_light - v(:,i);
    % Update electric and magnetic fields
    E(:,i) = (k*norm(R(:,i))./(dot(R(:,i),u(:,i)).^3)) .* ((speed_light^2 - norm(v(:,i))^2).*u(:,i) + cross( R(:,i), (cross(u(:,i),a(:,i))) ));
    B = 1/c^2 .* cross(v,E);
    % update acceleration and velocity to get next position of the charge
    acc(:,i) = q2.*(E(:,i) + cross(vel2(:,i),B(:,i)))/m2;
    vel2(:,i+1) = vel2(:,i) + acc(:,i).*dt;
    r_charge(:,i+1) = r_charge(:,i) + vel2(:,i).*dt;
end
% plot change in charge position
figure(6)
hold on
xlim([-4 4]);
ylim([-4 4]);
zlim([-8 8]);
plot3(r_charge(1,:),r_charge(2,:),r_charge(3,:));

