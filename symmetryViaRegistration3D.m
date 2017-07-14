function [planePoint, perpVector, scale] = symmetryViaRegistration3D(ptCloud)

% -------------------------
% pre-process

P = ptCloud;

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

% -------------------------
% reflect / register

VS = [1 0 0; 0 1 0; 0 0 1];
p = [0; 0; 0]; % point in plane
rmse = zeros(1,3);
for i = 1:3    
    V = VS(:,i);
    d = dot(p,V); % distance to origin
    S = [eye(3)-2*(V*V') 2*d*V; 0 0 0 1]; % symmetry transform
    Q = S*[P'; ones(1,size(P,1))];
    Q = Q(1:3,:)';
    fixed = pointCloud(P);
    moving = pointCloud(Q);
    [~,~,rmse(i)] = pcregrigid(moving, fixed, 'MaxIterations', 100);
end
[~,im] = min(rmse);
V = VS(:,im);
d = dot(p,V); % distance to origin
S = [eye(3)-2*(V*V') 2*d*V; 0 0 0 1]; % symmetry transform
N = V; % normal vector (used later)
Q = S*[P'; ones(1,size(P,1))];
Q = Q(1:3,:)';
fixed = pointCloud(P);
moving = pointCloud(Q);
tform = pcregrigid(moving, fixed, 'MaxIterations', 100);

% -------------------------
% compute symmetry plane

Tform = tform.T';
t = Tform(1:3,4);

T = S(1:3,1:3)*(Tform(1:3,1:3)');
[V,D] = eig(T);
ieig = [];
for i = 1:3
    if abs(D(i,i)+1) < 0.000001
        ieig = i;
        break;
    end
end
if isempty(ieig)
    error('no -1 eigenvalue')
end
V = V(:,ieig); % eigenvector of eigenvalue -1
p = (Tform(1:3,1:3)*(2*d*N)+t)/2; % point in plane

planePoint = p+m';
perpVector = V;

end