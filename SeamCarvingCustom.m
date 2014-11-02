function SeamCarvingCustom(FolderName , FileName , seamNumber , CROP, RATIO)
    warning('off','all');
    
    NKP = 10;

    Reader = ReadFunctions;
    Crop = CropFunctions;
    [video,nFrames,originalVidHeight,originalVidWidth,vidFPS,~,~] = Reader.ReadData(FolderName,FileName);
   
    mov = Reader.NewMovie(nFrames , originalVidHeight   ,originalVidWidth);
    mov = Reader.ReadMovie(mov , video );
    mov = Crop.RemoveBlackBars(mov);
    [vidHeight,vidWidth,~] = size(mov(1).cdata);

    % Define matrices for optimized seam carving
    seamVector = zeros(vidHeight , nFrames);
    KPE = zeros( NKP , nFrames);
    KPX = zeros( NKP , nFrames);
    KPY = zeros( NKP , nFrames);
    RPMAP = zeros(vidHeight , vidWidth);

    disp('Starting to remove seams...')
    for t=1:seamNumber; % t = number of seams removed
        
        % For the first frame, remove the seam without optimization
        EM = findEnergy(mov(1).cdata);
        seamImage = findSeamImg(EM);
        seamVector(:,1) = findSeam(seamImage);  
        mov(1).cdata = SeamCut(mov(1).cdata,seamVector(:,1));

        for n = 1:NKP

            start = floor((vidHeight/NKP)*(n-1))+1;
            eend  = floor((vidHeight/NKP)*(n));

            y = horzcat((start:eend)',seamVector(start:eend,1));
            energies = arrayfun( @( X )( EM( y( X , 1 ) , y( X , 2) ) ) , 1:length(y) );

            [val ind] = max( energies );
            KPE( n , 1 ) = val;
            KPX( n , 1 ) = y( ind , 2);
            KPY( n , 1 ) = y( ind , 1);

        end
        
        % For the remaining frames

        for k = 2 : nFrames

            EM = findEnergy(mov(k).cdata);
            seamImage = findSeamImg(EM);
            seamVector(:,k) = findSeam(seamImage);

            for n = 1:NKP

                start = floor((vidHeight/NKP)*(n-1))+1;
                eend  = floor((vidHeight/NKP)*(n));

                y = horzcat((start:eend)',seamVector(start:eend,1));
                energies = arrayfun( @( X )( EM( y( X , 1 ) , y( X , 2) ) ) , 1:length(y) );

                [val ind] = max( energies );
                KPE( n , k ) = val;
                KPX( n , k ) = y( ind , 2);
                KPY( n , k ) = y( ind , 1);

            end

            % calculate SAD for each KP (key point)
            RPMAP = zeros(vidHeight , vidWidth-t+1);
            for n = 1:NKP

                SR = vidHeight/NKP; % Search area constant
                MR = 3;             % Map area constant

                if( KPX(n,k-1)-MR < 1 )             MR = KPX(n,k-1)-1;                end
                if( KPX(n,k-1)+MR > vidWidth-t+1)   MR = vidWidth - KPX(n,k-1)-t+1;   end
                if( KPY(n,k-1)-MR < 1)              MR = KPY(n,k-1)-1;                end
                if( KPY(n,k-1)+MR > vidHeight)      MR = vidHeight - KPY(n,k-1);      end  

                if( KPX(n,k)-SR-MR < 1 )            SR = KPX(n,k)-MR-1;               end
                if( KPX(n,k)+SR+MR > vidWidth-t+1)  SR = vidWidth - KPX(n,k)-MR-t+1;  end
                if( KPY(n,k)-SR-MR < 1)             SR = KPY(n,k)-MR-1;               end
                if( KPY(n,k)+SR+MR > vidHeight)     SR = vidHeight - KPY(n,k)-MR;     end      

                if( SR < MR ) MR=SR;end
                if(SR==0 || MR == 0) continue; end
                % Define search area boundaries from previous KP's
                MA_X = KPX(n,k-1)-MR : KPX(n,k-1)+MR;
                MA_Y = KPY(n,k-1)-MR : KPY(n,k-1)+MR;

                SA_X = KPX(n,k)- SR : KPX(n,k)+SR;
                SA_Y = KPY(n,k)- SR : KPY(n,k)+SR;

                for x = 1:2*SR+1    

                    MA_SA_X = SA_X(x)-MR-x+1 : SA_X(x)+MR+x;
                    MA_SA_Y = SA_Y(x)-MR-x+1 : SA_Y(x)+MR+x;

                    RPMAP(MA_SA_Y( 1 ):MA_SA_Y( 2*MR+1 ),MA_SA_X( 1 ):MA_SA_X( 2*MR+1 )) = ...
                        imabsdiff( EM( MA_SA_Y( 1 ):MA_SA_Y( 2*MR+1 ),MA_SA_X( 1 ):MA_SA_X( 2*MR+1 ) ) ,...
                        EM( MA_Y( 1 ):MA_Y( 2*MR+1 ) , MA_X( 1 ):MA_X( 2*MR+1 ) ) );
                end
            end

            bot = 0;
            top = 0;
            top = floor(max(max(RPMAP)));
            bot = ceil(min(min(RPMAP(RPMAP>0.2))));
            if( ~(top==0 & bot==0) )
                ratio = 100/(top-bot);
                RPMAP = (RPMAP-bot) .* ratio;
            end
            RPMAP(RPMAP<0)=1;

            % Adjust EM and re-compute the seam
            EM = EM .* RPMAP;
            seamImage = findSeamImg(EM);
            seamVector(:,k) = findSeam(seamImage);
            mov(k).cdata = SeamCut(mov(k).cdata,seamVector(:,k));   
        end
        
        fprintf('\nNumber of seams removed = %d\n' , t);
    end

    for k = 1:nFrames
        mov(k).cdata = imresize(mov(k).cdata , [CROP.*RATIO(1,2) ,CROP.*RATIO(1,1)]);
    end
    movie2avi(mov, strcat(FolderName , FileName,'_seamCarved.avi') , 'fps',  vidFPS);

end

