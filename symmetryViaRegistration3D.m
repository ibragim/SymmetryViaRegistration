function [planePoints, perpVectors, strenghts, scale] = symmetryViaRegistration3D(P)
% [planePoints, perpVectors, strenghts, scale] = symmetryViaRegistration3D(ptCloud)
% Computes 3 guesses of symmetry planes on a 3D point cloud via registration.
%
% Input
%     P: a matrix NPx3 where NP is the number of points
%         and the 1st/2nd/3rd columns are their x/y/z coordinates.
%
% Outputs:
%     planePoints: a 3x3 matrix where planePoints(:,i) is a point in symmetry plane i
%     perpVectors: a 3x3 matrix where perpVectors(:,i) is a vector perpendicular to symmetry plane i
%     strenghts: a 1x3 vector where strenghts(i) is the strenght of sym. plane i;
%         symmetry plane parameters are output sorted in descending order of strength.
%     scale: equals max(std(ptCloud)); can be used to plot symmetry plane; see example below.
%
% Example
%     ptCloud = pcread('teapot.ply');
%     P = ptCloud.Location;
% 
%     [planePoints, perpVectors, strenghts, scale] = symmetryViaRegistration3D(P);
% 
%     p = planePoints(:,1);
%     v = perpVectors(:,1);
% 
%     % circle on symmetry plane centered at 'p' with radius 'scale'
%     w = [rand; rand; rand];
%     u = cross(v,w);
%     u = u/norm(u);
%     w = cross(v,u);
%     ags = 0:0.1:2*pi;
%     qs = zeros(3,length(ags));
%     for i = 1:length(ags)
%         ag = ags(i);
%         q = p+scale*(cos(ag)*u+sin(ag)*w);
%         qs(:,i) = q;
%     end
% 
%     figure
%     plot3(P(:,1),P(:,2),P(:,3),'r.','MarkerSize',1),hold on
%     plot3(p(1),p(2),p(3),'k*')
%     plot3([p(1) p(1)+scale/2*v(1)],[p(2) p(2)+scale/2*v(2)],[p(3) p(3)+scale/2*v(3)],'k-')
%     fill3(qs(1,:),qs(2,:),qs(3,:),'k'), alpha(0.1)
%     hold off
%     grid on, axis equal, axis vis3d
%     xlabel('x'), ylabel('y'), zlabel('z')
%     axis off
%     view(v)
%     ax = gca;
%     ax.Projection = 'perspective';
%     title('3D')
%     rotate3d on
%
% Reference
%     Finding Mirror Symmetry via Registration
%     Marcelo Cicconet, David G. C. Hildebrand, Hunter Elliott
%     https://arxiv.org/abs/1611.05971
%
% Marcelo Cicconet, Jul 2017


% pre-process

np = size(P,1);
if np > 10000
    rows = rand(1,np) < 10000/np;
    P = P(rows,:);
end

m = mean(P);
s = std(P);

scale = max(s);

P = P-repmat(m,[size(P,1) 1]);
P = P./repmat(s,[size(P,1) 1]);

% VS = [1 0 0; 0 1 0; 0 0 1];
VS = double(pca(P));
p = [0; 0; 0]; % point in plane
rmse = zeros(1,3);
planePoints = zeros(3,3);
perpVectors = zeros(3,3);
for i = 1:3    
    % reflect / register
    
    V = VS(:,i);
    d = dot(p,V); % distance to origin
    S = [eye(3)-2*(V*V') 2*d*V; 0 0 0 1]; % symmetry transform
    N = V; % normal vector (used later)
    Q = S*[P'; ones(1,size(P,1))];
    Q = Q(1:3,:)';
    fixed = pointCloud(P);
    moving = pointCloud(Q);
    [tform,~,rmse(i)] = pcregrigid(moving, fixed, 'MaxIterations', 100);

    % compute symmetry plane

    Tform = tform.T';
    t = Tform(1:3,4);

    T = S(1:3,1:3)*(Tform(1:3,1:3)');
    [V,D] = eig(T);
    jeig = [];
    for j = 1:3
        if abs(D(j,j)+1) < 0.000001
            jeig = j;
            break;
        end
    end
    if isempty(jeig)
        error('no -1 eigenvalue')
    end
    V = V(:,jeig); % eigenvector of eigenvalue -1
    p = (Tform(1:3,1:3)*(2*d*N)+t)/2; % point in plane

    planePoints(:,i) = p+m';
    perpVectors(:,i) = V;
end

strenghts = 1./rmse;
strenghts = strenghts/max(strenghts);

[strenghts, istrenghts] = sort(strenghts,'descend');
planePoints = planePoints(:,istrenghts);
perpVectors = perpVectors(:,istrenghts);

end