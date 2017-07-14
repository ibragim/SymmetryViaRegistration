clear, clc

scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
figure('Position',[scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2])

%% 2D

Image = imread('gantrycrane.png');
Image = imrotate(Image,45,'crop');

[angles, midPoints, segLengths] = symmetryViaRegistration2D(Image);

ag = angles(1);
mp = midPoints(:,1);
sl = segLengths(1);
p = mp+sl/2*[cos(ag); sin(ag)];
q = mp-sl/2*[cos(ag); sin(ag)];
Image = insertShape(Image,'line',[p(2) p(1) q(2) q(1)],'LineWidth',3,'Color','green');

subplot(1,2,1)
image(Image), axis equal, axis off, title('2D')


%% 3D

ptCloud = pcread('teapot.ply');
P = ptCloud.Location;

[planePoints, perpVectors, ~, scale] = symmetryViaRegistration3D(P);

p = planePoints(:,1);
v = perpVectors(:,1);

% circle on symmetry plane centered at 'p' with radius 'scale'
w = [rand; rand; rand];
u = cross(v,w);
u = u/norm(u);
w = cross(v,u);
ags = 0:0.1:2*pi;
qs = zeros(3,length(ags));
for i = 1:length(ags)
    ag = ags(i);
    q = p+scale*(cos(ag)*u+sin(ag)*w);
    qs(:,i) = q;
end

subplot(1,2,2)
plot3(P(:,1),P(:,2),P(:,3),'r.','MarkerSize',1),hold on
plot3(p(1),p(2),p(3),'k*')
plot3([p(1) p(1)+scale/2*v(1)],[p(2) p(2)+scale/2*v(2)],[p(3) p(3)+scale/2*v(3)],'k-')
fill3(qs(1,:),qs(2,:),qs(3,:),'k'), alpha(0.1)
hold off
grid on, axis equal, axis vis3d
xlabel('x'), ylabel('y'), zlabel('z')
axis off
view(v)
ax = gca;
ax.Projection = 'perspective';
title('3D')
rotate3d on