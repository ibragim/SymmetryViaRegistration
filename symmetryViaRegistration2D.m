function [angles, midPoints, segLengths, IOut] = symmetryViaRegistration2D(Image,varargin)

ip = inputParser;
ip.addParameter('RegMethod','dense'); % 'dense','sparse','nxc'
ip.addParameter('BoxSize',50); % only used with 'nxc' registration
ip.addParameter('NumBoxSamples',100); % only used with 'nxc' registration
ip.addParameter('MaxNumOutputs',1); % only used with 'nxc' registration
ip.parse(varargin{:});
prmts = ip.Results;

if nargout > 3
    IOut = Image;
end

I = Image;
if size(I,3) == 3
    I = rgb2gray(I);
elseif size(I,3) == 2 || size(I,3) > 3
    I = I(:,:,1);
end
ms = max(size(I));
newImSize = 200;
rf = newImSize/ms;
if rf < 1
    I = imresize(I,rf);
else
    rf = 1;
end

GI = normalize(imgradient(I));
I = normalize(double(I));

if strcmp(prmts.RegMethod,'nxc')
    [cellp,cellv] = eigsymNXC(GI,0,prmts.BoxSize,prmts.NumBoxSamples,prmts.MaxNumOutputs);
elseif strcmp(prmts.RegMethod,'dense') || strcmp(prmts.RegMethod,'sparse')
    [cellp,cellv] = eigsym(I,0,prmts.RegMethod);
else
    error('Registration method not recognized.')
end

nOutputs = length(cellp);

angles = zeros(1,nOutputs);
midPoints = zeros(nOutputs,2);
segLengths = zeros(1,nOutputs);

for iOutput = 1:nOutputs
    p = cellp{iOutput};
    v = cellv{iOutput};

    ag = atan2(v(1),v(2));
    if ag < 0
        ag = ag+pi;
    end

    xy = [];
    for j = 1:round(sqrt(2*newImSize^2))
        x = round(p(2)+j*cos(ag));
        y = round(p(1)+j*sin(ag));
        if x >= 1 && x <= size(I,1) && y >= 1 && y <= size(I,2)
            xy = [xy; [x y]];
        end
        x = round(p(2)-j*cos(ag));
        y = round(p(1)-j*sin(ag));
        if x >= 1 && x <= size(I,1) && y >= 1 && y <= size(I,2)
            xy = [xy; [x y]];
        end
    end
    if isempty(xy)
        xy = round(size(I)/2);
    else
        xy = round(mean(xy));
    end
    
    [midpointI,seglenI,~] = endpoints(I,1,ag,xy,'run');
    midpointI = midpointI/rf;
    seglenI = seglenI/rf;
    
    angles(iOutput) = ag;
    midPoints(iOutput,:) = midpointI;
    segLengths(iOutput) = seglenI;
    
    if nargout > 3
        v1 = midpointI+seglenI/2*[cos(ag) sin(ag)];
        v2 = midpointI-seglenI/2*[cos(ag) sin(ag)];
        IOut = insertShape(IOut,'line',[v1(2) v1(1) v2(2) v2(1)],'LineWidth',7,'Color','green');
        IOut = insertShape(IOut,'line',[v1(2) v1(1) v2(2) v2(1)],'LineWidth',5,'Color','yellow');
    end
end

end

% ----------------------------------------------------------------------------------------------------

function [p,v] = eigsymNXC(I,refAngle,boxSize,nBoxSamples,maxNOutputs)
% I should be double, in range [0,1]

% reflect
p = [size(I,2)/2; size(I,1)/2];
N = [cos(refAngle); sin(refAngle)];
d = dot(p,N);
S = [eye(2)-2*(N*N') 2*d*N; 0 0 1];

tform = affine2d(S');
xWorldLimits = [1 size(I,2)];
yWorldLimits = [1 size(I,1)];
J = imwarp(I,tform,'OutputView',imref2d(size(I),xWorldLimits,yWorldLimits));

% register
tforms = computeNormxcorrTransforms(J,I,boxSize,nBoxSamples,maxNOutputs);

nTForms = length(tforms);
p = cell(1,nTForms);
v = cell(1,nTForms);
for itform = 1:nTForms
    tform = tforms{itform};

    % compute sym line
    R = tform.T';
    t = R(1:2,3);

    T = S(1:2,1:2)*(R(1:2,1:2)');
    [V,D] = eig(T);
    ieig = [];
    for i = 1:2
        if abs(D(i,i)+1) < 0.000001
            ieig = i;
            break;
        end
    end
    if isempty(ieig)
        error('no -1 eigenvalue')
    end

    v{itform} = V(:,ieig); % eigenvector of eigenvalue -1
    v{itform} = [-v{itform}(2) v{itform}(1)]; % perp

    % point in line
    p{itform} = ((R(1:2,1:2)*(2*d*N))'+t')/2;
end

end

function tforms = computeNormxcorrTransforms(J,I,boxSize,nBoxSamples,maxNOutputs) % J: moving, I: fixed

[numrows,numcols] = size(I);
rangles = 0:6:360-6;
rmags = zeros(1,length(rangles));
vs = zeros(length(rangles),2);

parfor iangle = 1:length(rangles)
%     if mod(iangle,round(length(rangles)/10)) == 1
%         fprintf('.')
%     end
    A = zeros(2*numrows,2*numcols);
    flipI = imrotate(J,rangles(iangle),'crop');
    for index = 1:nBoxSamples
        w = boxSize;
        h = w;
        x0 = floor((size(flipI,2)-w)*rand);
        y0 = floor((size(flipI,1)-h)*rand);
        xcFlipI = x0+w/2;
        ycFlipI = y0+h/2;

        subFlipI = imcrop(flipI,[x0 y0 w h]);
        if var(subFlipI(:)) > 0.001
            [ROI,~,mc] = locateSubset(subFlipI,I);
            if mc > 0.25
                xcI = ROI(1)+ROI(3)/2;
                ycI = ROI(2)+ROI(4)/2;

                v = [xcI ycI]-[xcFlipI ycFlipI]; % translation

                row = max(min(round(v(2)+numrows),2*numrows),1);
                col = max(min(round(v(1)+numcols),2*numcols),1);
                A(row,col) = A(row,col)+1;
            end
        end
    end
    A = imfilter(A,fspecial('gaussian',12,3));
    maxA = max(A(:));
    
    [r,c] = find(A == maxA);
    v(1) = c(1)-numcols;
    v(2) = r(1)-numrows;
    rmags(iangle) = maxA;
    vs(iangle,:) = v;
end
% fprintf('\n')

[~,iangles] = sort(rmags,'descend');
nTForms = min(length(rangles),maxNOutputs);
tforms = cell(1,nTForms);
for i = 1:nTForms
    iangle = iangles(i);
    
    v = vs(iangle,:);

    arad = -rangles(iangle)/360*2*pi;

    % rotation with respect to center
    T1 = [eye(2) [-numcols/2; -numrows/2]; 0 0 1];
    T2 = [[cos(arad) -sin(arad); sin(arad) cos(arad)] [0; 0]; 0 0 1];
    T3 = [eye(2) [numcols/2; numrows/2]; 0 0 1];

    % translation
    T4 = [eye(2) v'; 0 0 1];

    % transform
    tform = affine2d((T4*T3*T2*T1)');
    
    tforms{i} = tform;
end

end

% ----------------------------------------------------------------------------------------------------

function [ROI,c,mc] = locateSubset(subI,I)

c = normxcorr2(subI,I);

mc = max(c(:));
[ypeak, xpeak] = find(c==mc);

yoffSet = ypeak(1)-size(subI,1);
xoffSet = xpeak(1)-size(subI,2);

ROI = [xoffSet+1, yoffSet+1, size(subI,2), size(subI,1)];

end

% ----------------------------------------------------------------------------------------------------

function [cellp,cellv] = eigsym(I,refAngle,regType)
% I should be double, in range [0,1]

% reflect
p = [size(I,2)/2; size(I,1)/2];
N = [cos(refAngle); sin(refAngle)];
d = dot(p,N);
S = [eye(2)-2*(N*N') 2*d*N; 0 0 1];

tform = affine2d(S');
xWorldLimits = [1 size(I,2)];
yWorldLimits = [1 size(I,1)];
J = imwarp(I,tform,'OutputView',imref2d(size(I),xWorldLimits,yWorldLimits));

if strcmp(regType,'dense')
%     optimizer = registration.optimizer.RegularStepGradientDescent;
%     optimizer.MaximumIterations = 500;
%     metric = registration.metric.MeanSquares;
    
    [optimizer,metric] = imregconfig('monomodal');
    tform = imregtform(J,I,'rigid',optimizer,metric);
elseif strcmp(regType,'sparse')
    metricThreshold = 10;
    numOctaves = 4;
    numScaleLevels = 4;

    pointsI = detectSURFFeatures(I,'MetricThreshold',metricThreshold,'NumOctaves',numOctaves,'NumScaleLevels',numScaleLevels);
    [featuresI, pointsI] = extractFeatures(I, pointsI, 'SURFSize', 128);

    pointsJ = detectSURFFeatures(J,'MetricThreshold',metricThreshold,'NumOctaves',numOctaves,'NumScaleLevels',numScaleLevels);
    [featuresJ, pointsJ] = extractFeatures(J, pointsJ, 'SURFSize', 128);

    indexPairs = matchFeatures(featuresJ, featuresI, 'Unique', true);%, 'MaxRatio', 0.1);

    matchedPointsI = pointsI(indexPairs(:,2), :);
    matchedPointsJ = pointsJ(indexPairs(:,1), :);

    if matchedPointsI.Count < 2 % not enough to apply estimateGeometricTransform
        tform = affine2d(eye(3,3));
    else
        tform = estimateGeometricTransform(matchedPointsJ, matchedPointsI,...
            'similarity', 'Confidence', 99.9, 'MaxNumTrials', 2000, 'MaxDistance', 10);
    end
    R = tform.T(1:2,1:2);
    R = 1/sqrt(det(R))*R; % get rid of 'scale' component
    tform.T(1:2,1:2) = R;
end

% compute sym line
R = tform.T';
t = R(1:2,3);

T = S(1:2,1:2)*(R(1:2,1:2)');
[V,D] = eig(T);
ieig = [];
for i = 1:2
    if abs(D(i,i)+1) < 0.000001
        ieig = i;
        break;
    end
end
if isempty(ieig)
    error('no -1 eigenvalue')
end

v = V(:,ieig); % eigenvector of eigenvalue -1
v = [-v(2) v(1)]; % perp

% point in line
p = ((R(1:2,1:2)*(2*d*N))'+t')/2;

cellp = {p};
cellv = {v};

end

% ----------------------------------------------------------------------------------------------------

function [midpointI,seglenI,proximity] = endpoints(I,sigma,angleI,midpointI,mode)

[nr,nc] = size(I);

imcent = round([nr nc]/2);
d = imcent-midpointI;
I = imtranslate(I,[d(2) d(1)]);
I = imrotate(I,-180*angleI/pi,'crop');

freq = 1/sigma;
index = 0;
anglerange = [-pi/3 -pi/6 0 pi/6 pi/3];
Convs = complex(zeros(nr,nc,length(anglerange)));
for langle = pi/2+anglerange
    index = index+1;
    [mr,mi] = cmorlet(sigma,freq,langle,0);
    J = conv2(I,mr+1i*mi,'same');
    Convs(:,:,index) = J;
end

hnc = floor(nc/2);
SS = zeros(nr,hnc);
for i = 1:index
   L = Convs(:,1:hnc,i);
   R = Convs(:,end-hnc+1:end,length(anglerange)+1-i);

   R = fliplr(R);
   S = abs(L.*conj(R));

   SS = max(SS,S);
end

proximity = sum(sum(SS)); % proximity between half images

SortS = normalize(sort(SS,2,'descend'));
SortS = SortS(:,1:10);

s = sum(SortS,2);

% ---------- find endpoints via 'local maxima' -----------------
l = 5;
k1 = [ones(1,l) -ones(1,l)];
s1 = conv(s,k1,'same');
if strcmp(mode,'debug')
    subplot(1,2,1)
    plot(s,'linewidth',2), hold on
    plot(s1,'linewidth',2)
end
[pks,lcs] = findpeaks(max(s1,0),'NPeaks',1,'MinPeakHeight',0.05*max(abs(s1)));
i0 = lcs;
if strcmp(mode,'debug')
    plot(i0,pks,'ok','linewidth',2)
end
[pks,lcs] = findpeaks(max(-flipud(s1),0),'NPeaks',1,'MinPeakHeight',0.05*max(abs(s1)));
i1 = length(s)-lcs+1;
if strcmp(mode,'debug')
    plot(i1,-pks,'ok','linewidth',2), hold off, axis off
end
if strcmp(mode,'debug')
    subplot(1,2,2)
    imshow([normalize(SS) SortS I])
    pause
end
% --------------------------------------------------------------


% ---------- find endpoints via 'fitting step function' --------
% ns = s/max(s);
% l = length(s);
% ijl = [];
% for i = 1:l-1
%     for j = i+1:l
%         L = ns(1:i);
%         M = ns(i:j);
%         R = ns(j:l);
%         loss = median(L)-0.1*median(M)+median(R);
%         ijl = [ijl; [i j loss]];
%     end
% end
% [~,im] = min(ijl(:,3));
% i0 = ijl(im,1);
% i1 = ijl(im,2);
% if strcmp(mode,'debug')
%     subplot(1,2,1)
%     plot(ns,'b'), hold on
%     x = 1:i0; plot(x,zeros(1,length(x)),'ok')
%     x = i0:i1; plot(x,ones(1,length(x)),'ok')
%     x = i1:l; plot(x,zeros(1,length(x)),'ok'), hold off
%     subplot(1,2,2)
%     imshow(imresize([normalize(SS) SortS I],2))
%     pause
% end
% --------------------------------------------------------------

% ---------- find endpoints via registration -------------------
% L = I(:,1:hnc);
% R = fliplr(I(:,end-hnc+1:end));
% s = 5;
% hh = 2;
% circBuffer = zeros(1,5);
% % i0 = nr-hh-1;
% for i = hh+1:nr-hh-1
%     LSlice = L(i-hh:i+hh,:);
%     RSlice = R(i-hh:i+hh,:);
%     count = 0;
%     for k = 1:hnc-s+1
%         j1 = k;
%         j2 = k+s-1;
%         Template = LSlice(:,j1:j2);
%         if var(Template(:)) > 0.001
%             [ROI,~,mc] = locate_subset(Template,RSlice);
%             if abs(ROI(1)-j1) < 10 && ROI(1) > 0 && ROI(1)+s-1 <= hnc && mc > 0% && max2(Template) > 0.1
%                 count = count+1;
%             end
%         end
%     end
%     circBuffer(mod(i-1,length(circBuffer))+1) = count;
%     if sum(circBuffer > 2) > 1
%         i0 = i;
%         break
%     end
% end
% % i1 = nr;
% for i = nr-hh:-1:i0
%     LSlice = L(i-hh:i+hh,:);
%     RSlice = R(i-hh:i+hh,:);
%     count = 0;
%     for k = 1:hnc-s+1
%         j1 = k;
%         j2 = k+s-1;
%         Template = LSlice(:,j1:j2);
%         if var(Template(:)) > 0.001
%             [ROI,~,mc] = locate_subset(Template,RSlice);
%             if abs(ROI(1)-j1) < 10 && ROI(1) > 0 && ROI(1)+s-1 <= hnc && mc > 0% && max2(Template) > 0.1
%                 count = count+1;
%             end
%         end
%     end
%     circBuffer(mod(i-1,length(circBuffer))+1) = count;
%     if sum(circBuffer > 2) > 1
%         i1 = i;
%         break
%     end
% end
% --------------------------------------------------------------

ep0 = [i0 round(nc/2)];
ep1 = [i1 round(nc/2)];

R = [cos(angleI) -sin(angleI); sin(angleI) cos(angleI)];
ep0 = round(R*(ep0-imcent)'+imcent'-d');
ep1 = round(R*(ep1-imcent)'+imcent'-d');

midpointI = round(0.5*(ep0+ep1))';
seglenI = norm(ep1-ep0);

end

% ----------------------------------------------------------------------------------------------------

function [mr,mi] = cmorlet(sigma,freq,angle,halfkernel)

% ref: http://arxiv.org/pdf/1203.1513.pdf, page 2

support = 2.5*sigma;

xmin = -support;
xmax = -xmin;
ymin = xmin;
ymax = xmax;
xdomain = xmin:xmax;
ydomain = ymin:ymax;

[x,y] = meshgrid(xdomain,ydomain);

xi = freq*[sin(angle); cos(angle)];

envelope = exp(-0.5*(x.*x+y.*y)/sigma^2); 
carrier = exp(1i*(xi(1)*x+xi(2)*y));

% makes sum of args = 0
C2 = sum(sum(envelope.*carrier))/sum(sum(envelope));

% makes sum of args*conj(args) = 1
arg = (carrier-C2).*envelope;
normfact = sum(sum(arg.*conj(arg)));
C1 = sqrt(1/normfact);

psi = C1*(carrier-C2).*envelope;

if halfkernel
    condition = ((xi(1)*x+xi(2)*y) <= 0);
    mr = real(psi).*condition;
    mi = imag(psi).*condition;
else
    mr = real(psi);
    mi = imag(psi);
end

end

% ----------------------------------------------------------------------------------------------------

function J = normalize(I)
    J = I-min(min(I));
    J = J/max(max(J));
end