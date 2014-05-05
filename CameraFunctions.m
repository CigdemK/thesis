function F = CameraFunctions
    F.CameraMotion = @CameraMotion;
end

function [H_out,e] = CameraMotion(mov)
    
    nFrames = length(mov);
    [fHeight, fWidth, ~] = size(mov(1).cdata);
    
    %---------------------------------------
    % Init functions (for jacobian)
    %---------------------------------------
    
    % Check if we already have the Jacobians for the given frame size
    varName = strcat('f_',num2str(fWidth),'_',num2str(fHeight));
    try
        functionJ = load('CameraMotionJacobians.mat',varName);
        foundFlag = 1;
    catch err
        foundFlag = 0;
    end
   
    if foundFlag == 0
        syms h11 h12 h13 h21 h22 h23 h31 h32 x y;
        varsVector = { h11; h12; h13; h21; h22; h23; h31; h32};
        
        F =  ((h11 * x + h12 * y+ h13 ) / ( h31 * x + h32 * y + 1 ) * ...
             ( h21 * x + h22 * y + h23 ) / ( h31 * x + h32 * y  + 1 ));
        functionJ = cell(8,1);
        for i = 1:8
            functionJ{i} =  matlabFunction(diff(F, varsVector{i} ) , 'vars' , varsVector);
        end
        
%         for h = 1: fHeight
%             H =  [H; [h12 * h , h32 * h , h22 * h , h32 * h ]];
%         end
%         
%         varsVector = { h11; h12; h13; h21; h22; h23; h31; h32};
%         for w = 1: fWidth
% %             F = @(h11, h12, h13, h21, h22, h23, h31, h32)...
% %                  ((h11 * w + h12 * h+ h13 ) / ( h31 * w + h32 * h + 1 ) * ...
% %                  ( h21 * w + h22 * h + h23 ) / ( h31 * w + h32 * h  + 1 ))
%             F = @(h11, h12, h13, h21, h22, h23, h31, h32)...
%                 ((h11 * w + H(:,1) + h13 ) ./ ( h31 * w + H(:,2) + 1 )) .* ...
%                  (( h21 * w + H(:,3) + h23 ) ./ ( h31 * w + H(:,4) + 1 ));
%             for i = 1:8
%                 functionJ{(w-1)*fHeight + h,i} =  diff(F, varsVector{i} );
%             end
%             toc;
%         end

%         varsVector = { h11; h12; h13; h21; h22; h23; h31; h32};

%         functionJ = cell(8,floor((fWidth*fHeight) / 20000)+1);
%         tic;
%         for i = 1:8
%             for k = 1:floor((fWidth*fHeight) / 20000) % the inner loop is for speeding up
%                 functionJ{i,k} =  diff( F((k-1)*20000+1:20000*k) , varsVector{i} );
%                 toc;
%             end
%             functionJ{i,31} =  diff( F(20000*k+1:end) , varsVector{i} );
%             
%         end
%         
%         functionJ = matlabFunction(functionJ , 'vars' , varsVector );
%         toc;
%         S.(varName) = functionJ;
%         if exist('CameraMotionJacobians.mat','file')
%             save('CameraMotionJacobians.mat', '-struct', 'S' , '-append');
%         else
%             save('CameraMotionJacobians.mat', '-struct', 'S');
%         end
%         toc;
    end
    
    %---------------------------------------
    % Calculate homography matrix for all frames
    %---------------------------------------
    Util = UtilFunctions;
    H = zeros(3,3,nFrames-1);
    H_out = zeros(9,nFrames-1);
    load('Panning0019.avi_H.mat','H','J');
    for t = 1:nFrames-1
%         H(:,:,t) = Util.Homography( mov(t).cdata, mov(t+1).cdata); % map frame t to frame t+1 
        H_tmp = H(:,:,t);
        H_out(:,t) = H_tmp(:);
    end
   
    %---------------------------------------
    % Adjust homogrpahy matrices iteratively
    %---------------------------------------
    ro = 1.1;
    nrIteration = 10;
    sigma = zeros(fHeight*fWidth,nFrames);
%     J = zeros(fHeight*fWidth,8,nFrames);
    tic;
    for k = 1:nrIteration 

        % Calculate sigma and J for all frames
        for t = 1:nFrames-1
            
            T = maketform( 'projective' , reshape( H_out(:,t) , [3 3] ) );
            img1 = mov(t).cdata;
            [imgT, ~, ~] = imtransform(img1, T , 'size' , size(img1));
            
            tmpSigma(:,:) = ((mov(t+1).cdata(:,:,1) - imgT(:,:,1)) + ...
                (mov(t+1).cdata(:,:,2) - imgT(:,:,2)) + ...
                (mov(t+1).cdata(:,:,3) - imgT(:,:,3))) / 3;
            sigma(:,t) = tmpSigma(:);

%             for w = 1: fWidth
%                 for h = 1: fHeight
%                     J((w-1)*fHeight + h,:,t) = arrayfun( @( X ) feval( functionJ{ X } ,... 
%                         H(1,1,t), H(1,2,t), H(1,3,t), ...   % variables for feval
%                         H(2,1,t), H(2,2,t), H(2,3,t), ...   
%                         H(3,1,t), H(3,2,t), w, h), ...      
%                         (1:8));                             % variables for arrayfun
%                 end
%             end

        end
        fprintf('Sigma & J calculation for iteration %d : %0.3f \n', k , toc);
        
        % Update dh and e iteratively
        dh = zeros(8,nFrames,nrIteration);
        lambda = zeros(1,nrIteration);
        mu = zeros(1,nrIteration)+0.1;
        e = sigma * -1; 
        
        for m = 1:nrIteration
            for t = 1:nFrames-1
                e(:,t+1) = S( mu(m) , ( J * dh(:,t,m) - sigma(:,t) + mu(m) * lambda(m) / 2 ) );
                dh(:,t,m+1) = J' / ( J' * J )' * ( sigma(:,t+1) + e(:,t+1) - lambda(m) / mu(m) );
                lambda(m+1) = lambda(m) + mu(m) * ( J * dh(:,t,m+1) - sigma(:,t+1) - e(:,t+1) );
            end
            mu(m+1) = ro*mu(m);
        end
        fprintf('dh & e calculation for iteration %d : %0.3f \n', k , toc);
        
        % Update H and e
        for t = 1:nFrames-1
            H_out(:,t) = H_out(:,t) + [ dh(:,t,m) ; 1 ];
        end
        
    end % endfor k

end

function [nr]= S(mu ,nr)

    e = nr-mu;
    e(e<0) = 0;
    
%     for i = 1:length(nr)
%         res = cell2mat(arrayfun(@(X)(mu*X + abs(X-e(i)) / 2) , 0:100, 'UniformOutput',false));
%         [~,ind] = min(res);
%         nr(i) = ind;
%     end
    nr=e;
    
end



