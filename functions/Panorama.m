function F = Panorama
    F.ImagePanorama = @ImagePanorama;
    F.MoviePanaroma = @MoviePanaroma;
end

function imgout = ImagePanorama(img1,img2)

    Motion = MotionFunctions;
    H = Motion.MyHomography(img1,img2);
   
    T = maketform('projective',H');
    [imgT, XDATA, YDATA] = imtransform(img2, T);
%     [imgT , XDATA, YDATA] = imtransform(img2, t, 'bicubic');

    % Stitch images
    [M1, N1, ~] = size(img1);
    [M2, N2, ~] = size(img2);
%     [MT, NT, ~] = size(imgT);

    % init masks
    mask1 = uint8(ones(size(img1, 1), size(img1, 2)));
    mask2 = uint8(ones(size(img2,1),size(img2,2)));
    mask2 = imtransform(mask2, T);   

    % stitched image bounds
    H=max( [M1 M1-YDATA(1) M2 M2+YDATA(1)] );
    W=max( [N1 N1-XDATA(1) N2 N2+XDATA(1)] );
    
    % Align image 1 bounds
    H1 = eye(3);
    if XDATA(1) < 0, H1(3,1)= -XDATA(1); end
    if YDATA(1) < 0, H1(3,2)= -YDATA(1); end
    T1 = maketform('affine',H1);

    [im1, ~, ~] = imtransform(img1, T1, 'XData', [1 W], 'YData', [1 H]);
    mask1 = imtransform(mask1, T1, 'XData', [1 W], 'YData', [1 H]);

    % Align image 2 bounds 
    H2 = eye(3);
    if XDATA(1) > 0, H2(3,1)= XDATA(1); end
    if YDATA(1) > 0, H2(3,2)= YDATA(1); end
    T2 = maketform('affine',H2);

    [im2, ~, ~] = imtransform(img2, T2, 'XData', [1 W], 'YData', [1 H]);
    mask2 = imtransform(mask2, T2, 'XData', [1 W], 'YData', [1 H]);

    % Size check
    if (size(im1,1) ~= size(im2,1)) || (size(im1,2) ~= size(im2,2))
        H = max( size(im1,1), size(im2,1) );
        W = max( size(im1,2), size(im2,2) );
        im1(H,W,:)=0;
        im2(H,W,:)=0;
        mask1(H,W)=0;
        mask2(H,W)=0;
    end

    % Combine both images
    n_layers = max(max(mask1));
    im1part = uint16(mask1 > (n_layers * mask2));
    im2part = uint16(mask2 > mask1);

    combpart = uint16(repmat(mask1 .* mask2,[1 1 3]));
    combmask = uint16(combpart > 0);

    stitched_image = repmat(im1part,[1 1 3]) .* uint16(im1) + repmat(im2part,[1 1 3]) .* uint16(im2);
    stitched_image = stitched_image + ( combpart .* uint16(im1) + combmask .* uint16(im2) ) ./ (combpart + uint16(ones(size(combpart,1),size(combpart,2),3)));
    imgout = uint8(stitched_image);
%     stitched_mask = mask1 + mask2;

%     % Blend images
%     Yoffset = 0;
%     Xoffset = 0;
%     pt = zeros(3,4);
%     [M1 N1 ~] = size(img1);
%     [M2 N2 ~] = size(img2);  
%     [M3 N3 ~] = size(imgT);
%     
%     % Find the location of the corner points of the transformed image
%     pt(:,1) = H*[1;1;1]; pt(:,2) = H*[N2;1;1];  pt(:,3) = H*[N2;M2;1]; pt(:,4) = H*[1;M2;1];
%     x2 = pt(1,:)./pt(3,:);
%     y2 = pt(2,:)./pt(3,:);
% 
%     up = round(min(y2));
%     left = round(min(x2));
%     if up  <=0; Yoffset = -up+1;    up = 1;   end
%     if left<=0; Xoffset = -left+1;  left = 1; end
% 
%     % Create panaroma image
%     
%     imgout(Yoffset+1:Yoffset+M1,Xoffset+1:Xoffset+N1,:) = img1;
%     imgout(up:up+M3-1,left:left+N3-1,:) = imgT;

end

function [panMov] = MoviePanaroma(moviePath )
    
    % Read data
    video=VideoReader(moviePath);
    nFrames = video.NumberOfFrames;
    frames = read(video);
    
    Cropper = CropFunctions;
    frames = Cropper.RemoveBlackBars(frames);
    
    % Create panaroma
    panMov = ImagePanorama(frames(:,:,:,1),frames(:,:,:,2));
    for n = 3:20
        panMov = ImagePanorama(panMov,mov(n).cdata);
    end
    
end
























