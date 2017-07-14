%% 2D

Image = imread('gantrycrane.png');
Image = imrotate(Image,45,'crop');

[angles, midPoints, segLengths, IOut] = symmetryViaRegistration2D(Image,'RegMethod','dense');

figure
image(IOut), axis equal, axis off, title('2D')

%% 3D

ptCloud = pcread('teapot.ply');
P = ptCloud.Location;

[planePoint, perpVector, scale] = symmetryViaRegistration3D(P);

p = planePoint;
v = perpVector;

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

figure
plot3(P(:,1),P(:,2),P(:,3),'r.','MarkerSize',1),hold on
plot3(p(1),p(2),p(3),'k*')
plot3([p(1) p(1)+scale/2*v(1)],[p(2) p(2)+scale/2*v(2)],[p(3) p(3)+scale/2*v(3)],'k-')
fill3(qs(1,:),qs(2,:),qs(3,:),'k'), alpha(0.1)
hold off
grid on, axis equal, axis vis3d
xlabel('x'), ylabel('y'), zlabel('z')
axis off
view(V)
ax = gca;
ax.Projection = 'perspective';
title('3D')


%% Reference

% Finding Mirror Symmetry via Registration
% Marcelo Cicconet, David G. C. Hildebrand, Hunter Elliott
% https://arxiv.org/abs/1611.05971