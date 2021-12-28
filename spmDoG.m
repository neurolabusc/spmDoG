function spmDoG(fnms, smallKernel, dehaze)
%Create a difference of Gaussian edge map
% fnms: filenames of images to process
% smallKernel: minimum kernel full-width half maximum (millimeters), default 4
% dehaze: zero dark voxels (1) or not (0), default 1
%
%Example 
% spmDoG('DWI.nii.gz') %use default FWHM and dehaze
% spmDoG('DWI.nii.gz', 2, 0) %2mm FWHM, no dehaze
% spmDoG; %use GUI
% spmDoG(strvcat('fMRI.nii','T1.nii'));

if ~exist('fnms','var')
	fnms = spm_select(inf,'^.*\.(gz|img|nii)$','Select image[s] for NaN removal'); 
    prompt = {'Small smooth FWHM (mm):','Dehaze (0=no, 1=yes)'};
    dlg_title = 'DoG options';
    def = {'4','1'};
    answer = inputdlg(prompt,dlg_title,1,def);
    smallKernel = str2num(answer{1});
    dehaze = str2num(answer{2});
end
if ~exist('smallKernel','var')
    smallKernel = 4;
end
if ~exist('dehaze','var')
    dehaze = 1;
end
bigKernel = smallKernel * 1.6;
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    fnm = nii_ungz(fnm);
    %load image
    hdr = spm_vol(fnm);
    img = spm_read_vols(hdr);
    % remove zero darkest voxels
    if (dehaze)
        img = removeHazeSub(img, 5);
    end
    %blur data: convert smoothing to mm, not voxels
    VOX = sqrt(sum(hdr.mat(1:3,1:3).^2));
    smallBlur  = smallKernel./VOX; 
    bigBlur  = bigKernel./VOX;
    img = differenceOfGaussian(img, smallBlur, bigBlur);
    %save results
    [pth,nam,ext, ~] = spm_fileparts(deblank(fnm));
    hfnm = fullfile(pth,['dog',  nam, ext]);
    hdr.fname = hfnm;
    hdr.pinfo(1) = 1/255;
    spm_write_vol(hdr,img);
end
end %matDoG()

function img = differenceOfGaussian(img, smallBlur, bigBlur)
    %small blur:
    simg = img + 0.0;
    spm_smooth(img,simg,smallBlur);
    %big blur
    bimg = img + 0.0; 
    spm_smooth(img,bimg,bigBlur);
    %difference of gaussians:
    img = simg - bimg;
    %binarize
    img = (img > 0) + 0;
    %outline
    %img = bwdist(img); %bwdist requires the image processing toolbox
    img = distance_field3D (img);

    img(img > 1) = 0;
    img(img > 0) = 1;
end % differenceOfGaussian()

function img = removeHazeSub(img, level)
%level = 1..5 (3/4,2/3,1/2,1/3,1/4
level = max(min(level,5),1);
classes = abs(level-3)+2;
ndark = max(4 - level, 1);
mask = otsuSub(img(:),classes);
mask = (mask > ndark) + 0;
fprintf('Otsu with %d classes, masking darkest %d, %g%% survive\n',classes, ndark, round(100.0 * sum(mask(:))/numel(mask),2) );
mask = reshape(mask, size(img));
img = img .* mask;
end %removeHazeSub()

function [IDX,sep] = otsuSub(I,n)
%OTSU Global image thresholding/segmentation using Otsu's method.
%   IDX = OTSU(I,N) segments the image I into N classes by means of Otsu's
%   N-thresholding method. OTSU returns an array IDX containing the cluster
%   indices (from 1 to N) of each point. Zero values are assigned to
%   non-finite (NaN or Inf) pixels.
%
%   IDX = OTSU(I) uses two classes (N=2, default value).
%
%   [IDX,sep] = OTSU(...) also returns the value (sep) of the separability
%   criterion within the range [0 1]. Zero is obtained only with data
%   having less than N values, whereas one (optimal value) is obtained only
%   with N-valued arrays.
%
%   Notes:
%   -----
%   It should be noticed that the thresholds generally become less credible
%   as the number of classes (N) to be separated increases (see Otsu's
%   paper for more details).
%
%   If I is an RGB image, a Karhunen-Loeve transform is first performed on
%   the three R,G,B channels. The segmentation is then carried out on the
%   image component that contains most of the energy.
%
%   Example:
%   -------
%   load clown
%   subplot(221)
%   X = ind2rgb(X,map);
%   imshow(X)
%   title('Original','FontWeight','bold')
%   for n = 2:4
%     IDX = otsu(X,n);
%     subplot(2,2,n)
%     imagesc(IDX), axis image off
%     title(['n = ' int2str(n)],'FontWeight','bold')
%   end
%   colormap(gray)
%
%   Reference:
%   ---------
%   Otsu N, <a href="matlab:web('http://dx.doi.org/doi:10.1109/TSMC.1979.4310076')">A Threshold Selection Method from Gray-Level Histograms</a>,
%   IEEE Trans. Syst. Man Cybern. 9:62-66;1979
%
%   See also GRAYTHRESH, IM2BW
%
%   -- Damien Garcia -- 2007/08, revised 2010/03
%   Visit my <a
%   href="matlab:web('http://www.biomecardio.com/matlab/otsu.html')">website</a> for more details about OTSU
% Copyright (c) 2010, Damien Garcia
error(nargchk(1,2,nargin))

% Checking n (number of classes)
if nargin==1
    n = 2;
elseif n==1;
    IDX = NaN(size(I));
    sep = 0;
    return
elseif n~=abs(round(n)) || n==0
    error('MATLAB:otsu:WrongNValue',...
        'n must be a strictly positive integer!')
elseif n>255
    n = 255;
    warning('MATLAB:otsu:TooHighN',...
        'n is too high. n value has been changed to 255.')
end

I = single(I);


% Convert to 256 levels
I = I-min(I(:));
I = round(I/max(I(:))*255);

% Probability distribution
unI = sort(unique(I));
nbins = min(length(unI),256);
if nbins==n
    IDX = ones(size(I));
    for i = 1:n, IDX(I==unI(i)) = i; end
    sep = 1;
    return
elseif nbins<n
    IDX = NaN(size(I));
    sep = 0;
    return
elseif nbins<256
    [histo,pixval] = hist(I(:),unI);
else
    [histo,pixval] = hist(I(:),256);
end
P = histo/sum(histo);
clear unI

% Zeroth- and first-order cumulative moments
w = cumsum(P);
mu = cumsum((1:nbins).*P);

% Maximal sigmaB^2 and Segmented image
if n==2
    sigma2B =...
        (mu(end)*w(2:end-1)-mu(2:end-1)).^2./w(2:end-1)./(1-w(2:end-1));
    [maxsig,k] = max(sigma2B);
    
    % segmented image
    IDX = ones(size(I));
    IDX(I>pixval(k+1)) = 2;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
elseif n==3
    w0 = w;
    w2 = fliplr(cumsum(fliplr(P)));
    [w0,w2] = ndgrid(w0,w2);
    
    mu0 = mu./w;
    mu2 = fliplr(cumsum(fliplr((1:nbins).*P))./cumsum(fliplr(P)));
    [mu0,mu2] = ndgrid(mu0,mu2);
    
    w1 = 1-w0-w2;
    w1(w1<=0) = NaN;
    
    sigma2B =...
        w0.*(mu0-mu(end)).^2 + w2.*(mu2-mu(end)).^2 +...
        (w0.*(mu0-mu(end)) + w2.*(mu2-mu(end))).^2./w1;
    sigma2B(isnan(sigma2B)) = 0; % zeroing if k1 >= k2
    
    [maxsig,k] = max(sigma2B(:));
    [k1,k2] = ind2sub([nbins nbins],k);
    
    % segmented image
    IDX = ones(size(I))*3;
    IDX(I<=pixval(k1)) = 1;
    IDX(I>pixval(k1) & I<=pixval(k2)) = 2;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
else
    k0 = linspace(0,1,n+1); k0 = k0(2:n);
    [k,y] = fminsearch(@sig_func,k0,optimset('TolX',1));
    k = round(k*(nbins-1)+1);
    
    % segmented image
    IDX = ones(size(I))*n;
    IDX(I<=pixval(k(1))) = 1;
    for i = 1:n-2
        IDX(I>pixval(k(i)) & I<=pixval(k(i+1))) = i+1;
    end
    
    % separability criterion
    sep = 1-y;
    
end

IDX(~isfinite(I)) = 0;


% Function to be minimized if n>=4
    function y = sig_func(k)
        
        muT = sum((1:nbins).*P);
        sigma2T = sum(((1:nbins)-muT).^2.*P);
        
        k = round(k*(nbins-1)+1);
        k = sort(k);
        if any(k<1 | k>nbins), y = 1; return, end
        
        k = [0 k nbins];
        sigma2B = 0;
        for j = 1:n
            wj = sum(P(k(j)+1:k(j+1)));
            if wj==0, y = 1; return, end
            muj = sum((k(j)+1:k(j+1)).*P(k(j)+1:k(j+1)))/wj;
            sigma2B = sigma2B + wj*(muj-muT)^2;
        end
        y = 1-sigma2B/sigma2T; % within the range [0 1]
        
    end
end %OtsuSub

function df = distance_field3D (img)
% compute distance field for a three-dimensional binary image.
% each voxel in the distance field is assigned to the distance to the
% nearest "true" voxel of the binary image.
% the Python code for two-dimensional image is at Philip Rideout's blog:
% https://prideout.net/blog/distance_fields/
% the pseudocode for one-dimensional vector is here:
% Felzenszwalb, P. F., & Huttenlocher, D. P. (2012). Distance transforms of sampled functions. Theory of computing, 8(1), 415-428.
% Matlab translation by Grigori Yourganov
img (img == 0) = Inf;
img (img == 1) = 0;
for i = 1:size (img, 2)
    for j = 1:size (img, 3)
        img1 (:, i, j) = horizontal_pass (img (:, i, j));
    end
end
img2 = permute (img1, [2 1 3]);
for i = 1:size (img2, 2)
    for j = 1:size (img2, 3)
        img3 (:, i, j) = horizontal_pass (img2 (:, i, j));
    end
end
img2 = permute (img3, [3 2 1]);
clear img3
for i = 1:size (img2, 2)
    for j = 1:size (img2, 3)
        img3 (:, i, j) = horizontal_pass (img2 (:, i, j));
    end
end
img3 = permute (img3, [2 3 1]);
df = sqrt (img3);
end % distance_field3D()

function df = horizontal_pass (single_row)
[hull_vertices, hull_intersections] = find_hull_parabolas (single_row);
df = march_parabolas (single_row, hull_vertices, hull_intersections);
end % horizontal_pass()

function [hull_vertices, hull_intersections] = find_hull_parabolas (single_row)
k = 1;
v(1) = 1;
z(1) = -Inf;
z(2) = Inf;
for q = 2:length (single_row)
    s = intersect_parabolas (single_row, q, v(k));
    while s <= z(k) && k > 1
        k = k - 1;
        s = intersect_parabolas (single_row, q, v(k));
    end
    k = k + 1;
    v(k) = q;
    z(k) = s;
    z(k+1) = Inf;
end
hull_vertices = v;
hull_intersections = z;
end %find_hull_parabolas()

function s = intersect_parabolas (single_row, p, q)
s = ((single_row (p) + p*p) - (single_row (q) + q*q)) / (2*p - 2*q);
if isnan (s)
    s = Inf;
end
end %intersect_parabolas()

function df = march_parabolas (single_row, hull_vertices, hull_intersections)
d = single_row;
v = hull_vertices;
z = hull_intersections;
k = 1;
for q = 1:length (d)
    while z(k+1) < q
        k = k + 1;
    end
    dx = q - v(k);
    single_row(q) = dx*dx + single_row(v(k));
end
df = single_row;
end % march_parabolas()

