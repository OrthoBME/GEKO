% function GEKO_v3(IDprefix,pics)

% GEKO is an automated GUI to evaluate knee OA in frontally sectioned
% rodent histological images. It calculates 8 of the 10 OARSI recommended
% grades from the OARSI histopathology initiative (Gerwin 2010).

IDprefix = 'mmtmag rat';
pics = [1 3 7];



lines = 'hide'; % 'hide' or 'show' to have lines appear as you draw them

% Import Reference Images
OARSI1 = imread('OARSI1.png');
OARSI2 = imread('OARSI2.png');
OARSI3 = imread('OARSI3.png');
OARSI4 = imread('OARSI4.png');
OARSI5 = imread('OARSI5.png');
OARSI6 = imread('OARSI6.png');
OARSI7 = imread('OARSI7.png');
OARSI8 = imread('OARSI8.png');
OARSI9 = imread('OARSI9.png');



for p = 1:length(pics)
%     clearvars( '-except', 'masterhist', 'p', 'pics', 'OARSI1', 'OARSI2', 'OARSI3', 'OARSI4', 'OARSI5', 'OARSI6', 'OARSI7', 'OARSI8', 'OARSI9', 'lines', 'IDprefix','wfile');
    
    tstart = fix(clock);
    
    display(['Processing "' IDprefix '_' int2str(pics(p)) '"'])
    
    
    %% FILE NAMING
% % %     % Automatically opens the image file if it's in the pathway. You can
% % %     % change line 37 (creation of variable "filename") to fit alternative
% % %     % file names. Also, letters can be used instead by substituting a
% % %     % string call in place of "int2str(pics(p))" in line 37.
%     filename = [ IDprefix int2str(pics(p)) '.jpg'];
    filename = [IDprefix int2str(pics(p)) 'R.jpg'];
%     if ~exist(filename)
%         filename(end-4:end) = '.png';
%         if ~exist(filename)
%             filename(end-4:end) = '.tiff';
%         end
%     end
%     
    if ~exist(filename)
        display('ERROR: Image file not found. \nPlease use and image with .jpg, .png, or .tiff extensions')
    end
    



    %write file name
    wfile = ['GEKO_' IDprefix '.xls'];
    
    
    I = imread(filename);
    ID = [IDprefix int2str(tstart(2)) '-' int2str(tstart(3)) '-' int2str(tstart(1)) '_IM_' int2str(pics(p))];
    
    [r c d] = size(I);
    
  
    
    %% Create the marking points
    redo = 'y';
    while redo == 'y';
        
        redo = 'n';
        
%         hFig = figure(1);
%         set(hFig, 'Position', [400 150 1000 700])
        figure(1)
        hold on
        subplot(3,4,1), image(OARSI1); axis off
        title('Reference Image')
        axis image
        subplot(3,4,[2:4,6:8,10:12]), image(I)
        title(ID)
        axis image
        
        linethick = 2;
        
        
        %tibial plateau
        % h = msgbox('Please click the two endpoints defining the tibial plateau.','Tibial Plateau','none','','','replace');
        % uiwait(h)
        clear px1 py1 px2 py2
        t1 = subplot(3,4,5, 'replace');
        text(0,0,sprintf('Please click the \ntwo endpoints defining \nthe tibial plateau.'), 'Parent', t1); axis off
        
        [px1 py1] = ginput(1);
        hold on
        subplot(3,4,[2:4,6:8,10:12]), plot(px1,py1,'ko');
        
        if exist('px1')
            [px2 py2] = ginput(1);
            subplot(3,4,[2:4,6:8,10:12]), plot(px2,py2,'ko');
        end
        
        if lines == 'show'
            subplot(3,4,[2:4,6:8,10:12]), plat = plot([px1 px2],[py1 py2],'k');
            set(plat,'LineWidth',linethick);
        end
        
        
        %Joint Capsule
        clear sx1 sy1 sx2 sy2
        subplot(3,4,1), image(OARSI2), axis off
        title('Reference Image')
        axis image
        % h = msgbox('Please click the two endpoints defining the joint capsule.','Joint Capsule','none','','','replace');
        % uiwait(h)
        
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('Please click the \ntwo endpoints defining \nthe joint capsule.'), 'Parent', t1); axis off
        
        [sx1 sy1] = ginput(1);
        subplot(3,4,[2:4,6:8,10:12]), plot(sx1,sy1,'go');
        
        if exist('sx1')
            [sx2 sy2] = ginput(1);
            subplot(3,4,[2:4,6:8,10:12]), plot(sx2,sy2,'go');
        end
        
        if lines == 'show'
            subplot(3,4,[2:4,6:8,10:12]), plat = plot([sx1 sx2],[sy1 sy2],'g');
            set(plat,'LineWidth',linethick);
        end
        
        
        %Growth Plates
        clear mgpx1 mgpy1 mgpx2 mgpy2 lgpx1 lgpy1 lgpx2 lgpy2
        subplot(3,4,1), image(OARSI3), axis off
        title('Reference Image')
        axis image
        % h = msgbox('Please click the two endpoints defining the thickness of the medial half of the medial compartment growth plate.','Growth Plate','None','','','replace');
        % uiwait(h)
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('Please click the \ntwo endpoints defining the \nthickness of the medial half \nof the medial compartment \ngrowth plate.'), 'Parent', t1); axis off
        
        [mgpx1 mgpy1] = ginput(1);
        subplot(3,4,[2:4,6:8,10:12]), plot(mgpx1,mgpy1,'mo');
        
        if exist('mgpx1')
            [mgpx2 mgpy2] = ginput(1);
            subplot(3,4,[2:4,6:8,10:12]), plot(mgpx2,mgpy2,'mo');
        end
        
        if lines == 'show'
            subplot(3,4,[2:4,6:8,10:12]), plat = plot([mgpx1 mgpx2],[mgpy1 mgpy2],'m');
            set(plat,'LineWidth',linethick);
        end
        
        
        
        subplot(3,4,1), image(OARSI4), axis off
        title('Reference Image')
        axis image
        % h = msgbox('Please click the two endpoints defining the thciness of the lateral half of the medial compartment growth plate.','Growth Plate','none','','','replace');
        % uiwait(h)
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('Please click the two endpoints defining the \nthickness of the lateral half \nof the medial compartment \ngrowth plate.'), 'Parent', t1); axis off
        
        [lgpx1 lgpy1] = ginput(1);
        subplot(3,4,[2:4,6:8,10:12]), plot(lgpx1,lgpy1,'mo');
        
        if exist('lgpx1')
            [lgpx2 lgpy2] = ginput(1);
            subplot(3,4,[2:4,6:8,10:12]), plot(lgpx2,lgpy2,'mo');
        end
        
        if lines == 'show'
            subplot(3,4,[2:4,6:8,10:12]), plat = plot([lgpx1 lgpx2],[lgpy1 lgpy2],'m:');
            set(plat,'LineWidth',linethick);
        end
        
        
        
        %Osteochondral Interface
        clear ocx1 ocy1 ocx2 ocy2
        subplot(3,4,1), image(OARSI5), axis off
        title('Reference Image')
        axis image
        % h = msgbox('Please trace the osteochondral interface and press enter when done.','Osteochondral Interface','none','','','replace');
        % uiwait(h)
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('Please trace the \nosteochondral interface and \npress enter when done.'), 'Parent', t1); axis off
        
        button = 1;
        i = 0;
        clear ocx ocy ocx1 ocy1
        while button == 1
            i = i+1;
            [ocx1 ocy1] = ginput(1);
            if size(ocx1,1)>0
                ocx(i) = ocx1;
                ocy(i) = ocy1;
                subplot(3,4,[2:4,6:8,10:12]), plot(ocx,ocy,'yo');
            else
                button = 0;
            end
        end
        
        if lines == 'show'
            subplot(3,4,[2:4,6:8,10:12]), oc = plot(ocx,ocy,'y');
            set(oc,'LineWidth',linethick);
        end
        
        
        %Total Cartilage Degeneration Width
        clear tcdwx1 tcdwy1 tcdwx2 tcdwy2
        subplot(3,4,1), image(OARSI6), axis off
        title('Reference Image')
        axis image
        % h = msgbox('If present, please click the two endpoints defining the total cartilage degeneration width. If not present, just press enter.','Total Cartilage Degeneration Width','none','','','replace');
        % uiwait(h)
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('If present, please click \nthe two endpoints defining \nthe total cartilage \ndegeneration width. \nIf not present, \njust press enter.'), 'Parent', t1); axis off
        
        [tcdwx1 tcdwy1] = ginput(1);
        if exist('tcdwx1')
            subplot(3,4,[2:4,6:8,10:12]), plot(tcdwx1,tcdwy1,'co');
        end
        
        if exist('tcdwx1')
            if size(tcdwx1)>0
                [tcdwx2 tcdwy2] = ginput(1);
                if exist('tcdwx1')
                    subplot(3,4,[2:4,6:8,10:12]), plot(tcdwx2,tcdwy2,'co');
                end
            end
        end
        
        if lines == 'show'
            if exist('tcdwx1')
                if size(tcdwx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot([tcdwx1 tcdwx2],[tcdwy1 tcdwy2],'c');
                    set(plat,'LineWidth',linethick);
                end
            end
        end
        
        %Significant Cartilage Degeneration Width
        clear scdwx1 scdwy1 scdwx2 scdwy2
        subplot(3,4,1), image(OARSI7), axis off
        title('Reference Image')
        axis image
        % h = msgbox('If present, please click the two endpoints defining the significant cartilage degeneration width. If not present, just press enter.','Significant Cartilage Degeneration Width','none','','','replace');
        % uiwait(h)
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('If present, please click \nthe two endpoints defining \nthe significant cartilage \ndegeneration width. \nIf not present, \njust press enter.'), 'Parent', t1); axis off
        
        [scdwx1 scdwy1] = ginput(1);
        if exist('scdwx1')
            subplot(3,4,[2:4,6:8,10:12]), plot(scdwx1,scdwy1,'bo');
        end
        
        if exist('scdwx1')
            if size(scdwx1)>0
                [scdwx2 scdwy2] = ginput(1);
                if exist('scdwx1')
                    subplot(3,4,[2:4,6:8,10:12]), plot(scdwx2,scdwy2,'bo');
                end
            end
        end
        
        if lines == 'show'
            if exist('scdwx1')
                if size(scdwx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot([scdwx1 scdwx2],[scdwy1 scdwy2],'b:');
                    set(plat,'LineWidth',linethick);
                end
            end
        end
        
        %Cartilage Lesion
        subplot(3,4,1), image(OARSI8), axis off
        title('Reference Image')
        axis image
        clear cx cy cx1 cy1
        % h = msgbox('If present, trace the cartilage lesion and press enter when done. If not present, just press enter.','Cartilage Lesion','none','','','replace');
        % uiwait(h)
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('If present, trace \nthe cartilage lesion \nand press enter when done. \nIf not present, \njust press enter.'), 'Parent', t1); axis off
        
        button = 1;
        i = 0;
        while button == 1
            i = i+1;
            [cx1 cy1] = ginput(1);
            if size(cx1,1)>0
                cx(i) = cx1;
                cy(i) = cy1;
                subplot(3,4,[2:4,6:8,10:12]), plot(cx,cy,'ro');
            else
                button = 0;
            end
        end
        
        if lines == 'show'
            if exist('cx')
                subplot(3,4,[2:4,6:8,10:12]), cart = fill(cx,cy,'r');
            end
        end
        
        
        %Osteophyte
        clear ostx1 osty1 ostx2 osty2
        subplot(3,4,1), image(OARSI9), axis off
        title('Reference Image')
        axis image
        % h = msgbox('If present, please click the two endpoints defining the thickest diameter of the osteophyte from the chondro-osseous junction to the cartilage surface. If not present, just press enter.','Osteophyte','none','','','replace');
        % uiwait(h)
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('If present, please click \nthe two endpoints defining \nthe thickest diameter \nof the osteophyte \nfrom the chondro-osseous junction \nto the cartilage surface. \nIf not present, \njust press enter.'), 'Parent', t1); axis off
        
        [ostx1 osty1] = ginput(1);
        if exist('ostx1')
            if size(ostx1)>0
                subplot(3,4,[2:4,6:8,10:12]), plot(ostx1,osty1,'yo');
            end
        end
        
        if exist('ostx1')
            if size(ostx1)>0
                [ostx2 osty2] = ginput(1);
                if exist('ostx1')
                    subplot(3,4,[2:4,6:8,10:12]), plot(ostx2,osty2,'yo');
                end
            end
        end
        
        if lines == 'show'
            if exist('ostx1')
                subplot(3,4,[2:4,6:8,10:12]), plat = plot([ostx1 ostx2],[osty1 osty2],'y:');
                set(plat,'LineWidth',linethick);
            end
        end
        
        
        
        %plotting lines if the use wants hidden lines
        if lines == 'hide'
            hold on
            
            subplot(3,4,[2:4,6:8,10:12]), plat = plot([px1 px2],[py1 py2],'k');
            set(plat,'LineWidth',linethick);
            
            
            subplot(3,4,[2:4,6:8,10:12]), plat = plot([sx1 sx2],[sy1 sy2],'g');
            set(plat,'LineWidth',linethick);
            
            
            subplot(3,4,[2:4,6:8,10:12]), plat = plot([mgpx1 mgpx2],[mgpy1 mgpy2],'m');
            set(plat,'LineWidth',linethick);
            
            
            subplot(3,4,[2:4,6:8,10:12]), plat = plot([lgpx1 lgpx2],[lgpy1 lgpy2],'m:');
            set(plat,'LineWidth',linethick);
            
            
            subplot(3,4,[2:4,6:8,10:12]), oc = plot(ocx,ocy,'y');
            set(oc,'LineWidth',linethick);
            
            
            if exist('tcdwx1')
                if size(tcdwx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot([tcdwx1 tcdwx2],[tcdwy1 tcdwy2],'c');
                    set(plat,'LineWidth',linethick);
                end
            end
            
            
            if exist('scdwx1')
                if size(scdwx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot([scdwx1 scdwx2],[scdwy1 scdwy2],'b:');
                    set(plat,'LineWidth',linethick);
                end
            end
            
            
            if exist('cx')
                if size(cx)>0
                    subplot(3,4,[2:4,6:8,10:12]), cart = fill(cx,cy,'r');
                end
            end
            
            
            if exist('ostx1')
                if size(ostx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot([ostx1 ostx2],[osty1 osty2],'y:');
                    set(plat,'LineWidth',linethick);
                end
            end
            
        end
        
        
        
        t1 = subplot(3,4,5,'replace');
        text(0,0,sprintf('If you would like to regrade the knee, \nplease press "y" in the command line, \nthen press enter. \n \nOtherwise, press enter.'), 'Parent', t1); axis off
        
        yesno = input('If you would like to regrade the knee, please press "y" in the command line, then press enter. Otherwise, press enter.', 's');
        if length(yesno)==0
            yesno = 'n';
        end
        redo = yesno;
        
    end
    
    
    
    % h = msgbox('Processing data. Please wait.','Processing','none','','','replace');
    t1 = subplot(3,4,5,'replace');
    text(0,0,sprintf('Processing data. \nPlease wait.'), 'Parent', t1); axis off
    
    
    saveas(figure(1), [ID '_GEKO_Figure.png']);
    
    
    
    
    
    
    
    
    %% 1. Cartilage Matrix Loss Width (0%, 50%, and 100% lesion depth)
    
    plateau = sqrt((px1-px2)^2+(py1-py2)^2);
    theta = atand((py1-py2)/(px1-px2));
    
    if exist('cx')
        
        %Determining whether tibial line (black) is above or below lesion)
        m = (py2-py1)/(px2-px1);
        b = py1-m*px1;
        lescent_x = mean(cx);
        lescent_y = mean(cy);
        newb = lescent_y-m*lescent_x;
        % comparison is reversed b/c image origin is in a different place than
        % cartesian plane.
        tibline = sign(newb-b);
        
        if cx(end-1)-cx(1)>cx(end)-cx(1)
            cx = cx(1:end-1);
            cy = cy(1:end-1);
        end
        
        
        h1 = sqrt((px1-cx(1))^2+(py1-cy(1))^2);
        h2 = sqrt((px2-cx(1))^2+(py2-cy(1))^2);
        h3 = sqrt((px1-px2)^2+(py1-py2)^2);
        
        s = (h1+h2+h3)/2;
        H = sqrt((4*s*(s-h1)*(s-h2)*(s-h3))/h3^2);
        
        dely = tibline*H*sind(90-theta);
        delx = tibline*H*cosd(90-theta);
        
        morey = tibline*round(r/30)*sind(90-theta);
        morex = tibline*round(r/30)*cosd(90-theta);
        
        
        oldplat = [px1 px2 py1 py2];
        
        px1 = (px1-delx)+(tibline)*morex;
        px2 = (px2-delx)+(tibline)*morex;
        py1 = (py1+dely)+(tibline*-1)*morey;
        py2 = (py2+dely)+(tibline*-1)*morey;
        
        
        
        
        % Creating the lesion on the image
        Ltemp = zeros(r,c);
        for i = 2:size(cx,2)
            %         rpts = linspace(r-cy(i-1),r-cy(i),1000);   %# A set of row points for the line
            rpts = linspace(cy(i-1),cy(i),1000);   %# A set of row points for the line
            cpts = linspace(cx(i-1),cx(i),1000);   %# A set of column points for the line
            index = sub2ind([r c],round(rpts),round(cpts));
            Ltemp(index) = 1;
            if i == size(cx,2)
                %             rpts = linspace(r-cy(1),r-cy(i),1000);   %# A set of row points for the line
                rpts = linspace(cy(1),cy(i),1000);   %# A set of row points for the line
                cpts = linspace(cx(1),cx(i),1000);   %# A set of column points for the line
                index = sub2ind([r c],round(rpts),round(cpts));
                Ltemp(index) = 1;
            end
        end
        
        for i = 2:size(ocx,2)
            %         rpts = linspace(r-ocy(i-1),r-ocy(i),1000);   %# A set of row points for the line
            rpts = linspace(ocy(i-1),ocy(i),1000);   %# A set of row points for the line
            cpts = linspace(ocx(i-1),ocx(i),1000);   %# A set of column points for the line
            index = sub2ind([r c],round(rpts),round(cpts));
            Ltemp(index) = 1;
        end
        
        stats = regionprops(Ltemp,'FilledImage','BoundingBox','Centroid','MajorAxisLength','Extrema', 'ConvexArea', 'ConvexHull','ConvexImage','Image','Orientation');
        lesion = imrotate(stats.FilledImage,theta);
        
        lesion = imdilate( lesion,strel('square',3));
        lesion =imerode(lesion,strel('square',2));
        [lesion lnum] = bwlabel(lesion(:,:,1));
        
        
        if lnum>1
            index = [min(find(sum(lesion==1,1))) min(find(sum(lesion==1,2)))];
            target = find(index==max(index));
            lesion = lesion(:,min(find(sum(lesion==target,1))):max(find(sum(lesion==target,1))));
            
            [csr csc] = find(lesion==target);
            cart_surf = min(csr);
            switch target
                case 1
                    [tmr tmc] = find(lesion==2);
                    tidemark = min(tmr);
                case 2
                    [tmr tmc] = find(lesion==1);
                    tidemark = min(tmr);
            end
            lesion = lesion(cart_surf:tidemark,:);
            lesion = lesion==target;
            
        else
            Ltemp = zeros(r,c);
            for i = 2:size(cx,2)
                rpts = linspace(cy(i-1),cy(i),1000);   %# A set of row points for the line
                cpts = linspace(cx(i-1),cx(i),1000);   %# A set of column points for the line
                index = sub2ind([r c],round(rpts),round(cpts));
                Ltemp(index) = 1;
                if i == size(cx,2)
                    rpts = linspace(cy(1),cy(i),1000);   %# A set of row points for the line
                    cpts = linspace(cx(1),cx(i),1000);   %# A set of column points for the line
                    index = sub2ind([r c],round(rpts),round(cpts));
                    Ltemp(index) = 1;
                end
            end
            stats = regionprops(Ltemp,'FilledImage','BoundingBox','Centroid','MajorAxisLength','Extrema', 'ConvexArea', 'ConvexHull','ConvexImage','Image','Orientation');
            lesion = imrotate(stats.FilledImage,theta);
            
            lesion = imdilate( lesion,strel('square',3));
            lesion =imerode(lesion,strel('square',2));
            [lesion lnum] = bwlabel(lesion(:,:,1));
            lesion = lesion(min(find(sum(lesion,2)>0)):max(find(sum(lesion,2)>0)),:);
            
        end
        
        
        zstats = regionprops(lesion,'FilledImage','BoundingBox','Centroid','MajorAxisLength','Extrema', 'ConvexArea', 'ConvexHull','ConvexImage','Image','Orientation');
        
        convexles = zstats.ConvexImage;
        lesstats = regionprops(convexles, 'Extrema');
        extrema = lesstats.Extrema;
        
        perc10 = round(.1*size(lesion,2));
        perc90 = round(.9*size(lesion,2));
        
        left = extrema(extrema(:,1)<perc10,:);
        left = left(left(:,2)==(min(left(:,2))),:);
        if size(left,1)>1
            left = left(left(:,1)==(min(left(:,1))),:);
            if size(left,1)>1
                left = left(1,:);
            end
        end
        
        right = extrema(extrema(:,1)>perc90,:);
        right = right(right(:,2)==(min(right(:,2))),:);
        if size(right,1)>1
            right = right(right(:,1)==(min(right(:,1))),:);
            if size(right,1)>1
                right = right(1,:);
            end
        end
        
        
        extrapatch = left(2)-min(zstats.Extrema(:,2));
        
        
        alpha = atand((left(2)-right(2))/(left(1)-right(1)));
        lesion = imrotate(lesion,.8*alpha);
        %     lesion = lesion(min(find(sum(lesion,2)>0)):max(find(sum(lesion,2)>0)),:);
        
        lesionbrackets = linspace(1,size(lesion,1),12);
        
        
        
        lesion0 = lesion(1:max([round(lesionbrackets(2)) round(left(2))+1]),:);
        lesion50 = lesion(round(lesionbrackets(6)):round(lesionbrackets(7)),:);
        lesion100 = lesion(round(lesionbrackets(11)):end,:);
        
            figure
            hold on
            subplot(3,2,3), imshow(lesion), title([ ID ' Cartilage Lesion'])
            hLine = imline(gca,[1 size(lesion,2)],[max([round(lesionbrackets(2)) round(left(2))+1]) max([round(lesionbrackets(2)) round(left(2))+1])]);
            hLine = imline(gca,[1 size(lesion,2)],[round(lesionbrackets(6)) round(lesionbrackets(6))]);
            hLine = imline(gca,[1 size(lesion,2)],[round(lesionbrackets(7)) round(lesionbrackets(7))]);
            hLine = imline(gca,[1 size(lesion,2)],[round(lesionbrackets(11)) round(lesionbrackets(11))]);
        
            subplot(3,2,2), imshow(lesion0), title('Cartilage Lesion at 0% Depth (Cartilage Surface)')
            subplot(3,2,4), imshow(lesion50), title('Cartilage Lesion at 50% Depth (Cartilage Midzone)')
            subplot(3,2,6), imshow(lesion100), title('Cartilage Lesion at 100% Depth (Tidemark)')
        
        %two approaches via solid rows and absolute columns:
        %widest row of lesion
        width0 = max(sum(lesion0,2));
        width50 = max(sum(lesion50,2));
        width100 = max(sum(lesion100,2));
        
        % %     %left and right most values of lesion within the image
        % %     width0_ = max(find((sum(lesion0)>0))) - min(find((sum(lesion0)>0)));
        % %     width50_ = max(find((sum(lesion50)>0))) - min(find((sum(lesion50)>0)));
        % %     width100_ = max(find((sum(lesion100)>0))) - min(find((sum(lesion100)>0)));
        
        
        
        
        
        
        %% Cartilage Degeneration Score
        lesfill = stats.FilledImage;
        Ltemp = zeros(r,c);
        Ltemp(round(stats.BoundingBox(2)):round(stats.BoundingBox(2))+size(lesfill,1)-1,round(stats.BoundingBox(1)):round(stats.BoundingBox(1))+size(lesfill,2)-1)=lesfill;
        %     Ltemp = Ltemp<1;
        justles = Ltemp;
        
        %     rpts = linspace(r-py1,r-py2,1000);   %# A set of row points for the line
        rpts = linspace(py1,py2,1000);
        cpts = linspace(px1,px2,1000);   %# A set of column points for the line
        index = sub2ind([r c],round(rpts),round(cpts));
        Ltemp(index) = 1;
        
        L = Ltemp;
        L = imdilate( L,strel('square',10));
        L =imerode(L,strel('square',8));
        [L lnum] = bwlabel(L(:,:,1));
        
        %     L = imrotate(L,theta+alpha);
        L = imrotate(L,theta);
        L = L(min(find(sum(L,2)>0)):max(find(sum(L,2)>0)),min(find(sum(L,1)>0)):max(find(sum(L,1)>0)));
        for i = 1:lnum
            holdmat=ones(size(L));
            holdmat(L(:)~=i)=0;
            
            cstats = regionprops(holdmat,'Centroid','Area','BoundingBox');
            Lmat(1,i) = i;
            Lmat(2,i) = cstats.Area;
            Lmat(3,i) = cstats.Centroid(1);
            Lmat(4,i) = cstats.Centroid(2);
            Lmat(5,i) = cstats.BoundingBox(1);
            Lmat(6,i) = cstats.BoundingBox(2);
            Lmat(7,i) = cstats.BoundingBox(3);
            Lmat(8,i) = cstats.BoundingBox(4);
        end
        cartles = Lmat(:,Lmat(7,:)==min(Lmat(7,:)));
        L = L==cartles(1);
        L = imrotate(L,alpha);
        L = L((min(find(sum(L,2)>0))):(max(find(sum(L,2)>0))),:);
        
        
        zones = linspace(1,size(L,2),4);
        
        zone1 = L(:,round(zones(3)):end);
        zone2 = L(:,round(zones(2)):round(zones(3)));
        zone3 = L(:,1:round(zones(2)));
        
        CDS(1) = max(sum(zone1,2))/round(zones(2)); ...zone1
            CDS(2) = max(sum(zone2,2))/round(zones(2)); ...zone2
            CDS(3) = max(sum(zone3,2))/round(zones(2)); ...zone3
            
        for i = 1:3
            if CDS(i)<0.05
                CDS(i) = 0;
            elseif CDS(i)<0.1
                CDS(i) = 1;
            elseif CDS(i)<0.25
                CDS(i) = 2;
            elseif CDS(i)<0.5
                CDS(i) = 3;
            elseif CDS(i)<0.75
                CDS(i) = 4;
            else
                CDS(i) = 5;
            end
        end
        CDSz1 = CDS(1);
        CDSz2 = CDS(2);
        CDSz3 = CDS(3);
        
        
        %%
            figure
            subplot(2,3,[1:3]), imshow(L), title([ ID ' Medial Tibial Compartment']);
            subplot(2,3,4), imshow(zone3), title('Zone 3 (Lateral Medial Compartment)');
            subplot(2,3,5), imshow(zone2), title('Zone 2 (Central Medial Compartment)');
            subplot(2,3,6), imshow(zone1), title('Zone 1 (Medial Medial Compartment)');
        %%
        
        
        
        
        
        %% Zonal Depth Ratio
        
        Ltemp(round(stats.BoundingBox(2)):round(stats.BoundingBox(2))+size(lesfill,1)-1,round(stats.BoundingBox(1)):round(stats.BoundingBox(1))+size(lesfill,2)-1)=lesfill;
        
        %     rpts = linspace(r-py1,r-py2,1000);   %# A set of row points for the line
        rpts = linspace(py1,py2,1000);   %# A set of row points for the line
        cpts = linspace(px1,px2,1000);   %# A set of column points for the line
        index = sub2ind([r c],round(rpts),round(cpts));
        Ltemp(index) = 1;
        
        for i = 2:size(ocx,2)
            %         rpts = linspace(r-ocy(i-1),r-ocy(i),1000);   %# A set of row points for the line
            rpts = linspace(ocy(i-1),ocy(i),1000);   %# A set of row points for the line
            cpts = linspace(ocx(i-1),ocx(i),1000);   %# A set of column points for the line
            index = sub2ind([r c],round(rpts),round(cpts));
            Ltemp(index) = 1;
        end
        
        L = Ltemp;
        L = imdilate( L,strel('square',8));
        L =imerode(L,strel('square',6));
        [L lnum] = bwlabel(L(:,:,1));
        
        
        L = imrotate(L,theta);
        L = L(min(find(sum(L,2)>0)):max(find(sum(L,2)>0)),min(find(sum(L,1)>0)):max(find(sum(L,1)>0)));
        L = L(6:end,:);
        L = L(min(find(sum(L,2)>0)):max(find(sum(L,2)>0)),min(find(sum(L,1)>0)):max(find(sum(L,1)>0)));
        
        
        [L lnum] = bwlabel(L);
        zones = linspace(1,size(L,2),7);
        
        zone1 = L(:,round(zones(5)):end);
        zone2 = L(:,round(zones(3)):round(zones(5)));
        zone3 = L(:,1:round(zones(3)));
        
        center = [zone1(:,round(zones(2))) zone2(:,round(zones(2))) zone3(:,round(zones(2)))];
        
        if lnum==1
            for i = 1:3
                if sum(center(:,i)==1)>0  % red=2, yellow=1
                    [holdmat hnum] = bwlabel(center(:,i));
                    if hnum>1
                        redmax = max(find(holdmat==1));
                        yellowmin = min(find(holdmat==2));
                        ZDR(1,i) = redmax;
                        ZDR(2,i) = yellowmin;
                        ZDR(3,i) = redmax/yellowmin;
                    else
                        ZDR(1,i) = 0;
                        ZDR(2,i) = 0;
                        ZDR(3,i) = 0;
                    end
                else
                    ZDR(1,i) = 0;
                    ZDR(2,i) = 0;
                    ZDR(3,i) = 0;
                end
            end
            ZDRz1 = ZDR(end,1);
            ZDRz2 = ZDR(end,2);
            ZDRz3 = ZDR(end,3);
            
            
        else
            for i = 1:3
                if sum(center(:,i)==1)>0 && sum(center(:,i)==2)>0 % red=2, yellow=1
                    redmax = max(find(center(:,i)==2));
                    yellowmin = min(find(center(:,i)==1));
                    ZDR(1,i) = redmax;
                    ZDR(2,i) = yellowmin;
                    ZDR(3,i) = redmax/yellowmin;
                else
                    ZDR(1,i) = 0;
                    ZDR(2,i) = 0;
                    ZDR(3,i) = 0;
                end
            end
            
            
            ZDRz1 = ZDR(end,1);
            ZDRz2 = ZDR(end,2);
            ZDRz3 = ZDR(end,3);
        end
        
        %%
        
            figure
        %     image(25*flipdim(L,1))
            image(25*L)
            hold on
            plot(zones(3)*ones(1,2),[1 size(L,1)],'r-')
            plot(zones(5)*ones(1,2),[1 size(L,1)],'r-')
            plot(zones(2)*ones(1,2),[1 size(L,1)],'y-.')
            plot(zones(4)*ones(1,2),[1 size(L,1)],'y-.')
            plot(zones(6)*ones(1,2),[1 size(L,1)],'y-.')
            title([ ID ' Medial Tibial Compartment']);
            legend('Zone Boundary','Zone Boundary','Zone1 Zonal Ratio Axis','Zone 2 Zonal Ratio Axis','Zone 3 Zonal Ratio Axis', 'Location','NorthWest')
        
            figure
            hold on
            image(25*flipdim(zone3,1))
            plot(zones(2)*ones(1,2),[1 size(zone3,1)],'y-.')
            title([ ID ' Zone 3 Medial Tibial Compartment']);
        
        
            figure
            hold on
            image(25*flipdim(zone2,1))
            plot(zones(2)*ones(1,2),[1 size(zone2,1)],'y-.')
            title([ ID ' Zone 2 Medial Tibial Compartment']);
        
            figure
            hold on
            image(25*flipdim(zone1,1))
            plot(zones(2)*ones(1,2),[1 size(zone1,1)],'y-.')
            title([ ID ' Zone 1 Medial Tibial Compartment']);
        %%
        
        
    else
        width0 = 0;
        width50 = 0;
        width100 = 0;
        
        CDSz1 = 0;
        CDSz2 = 0;
        CDSz3 = 0;
        
        ZDRz1 = 0;
        ZDRz2 = 0;
        ZDRz3 = 0;
        
    end
    
    
    %%   Total and  Significant and total Cartilage Degeneration Width
    
    if exist('tcdwx1')
        if size(tcdwx1)>0
            TCDW = sqrt((tcdwx1-tcdwx2)^2+(tcdwy1-tcdwy2)^2);
        else
            TCDW = 0;
        end
    else
        TCDW = 0;
    end
    if exist('scdwx1')
        if size(scdwx1)>0
            SCDW = sqrt((scdwx1-scdwx2)^2+(scdwy1-scdwy2)^2);
        else
            SCDW = 0;
        end
    else
        SCDW = 0;
    end
    
    %% 6. Osteophyte
    
    
    if exist('ostx1')
        if size(ostx1)>0
            OST = sqrt((ostx1-ostx2)^2+(osty1-osty2)^2);
        else
            OST = 0;
        end
    else
        OST = 0;
    end
    
    
    
    %% 7. Calcified Cartilage and Subchondral bone Damage Score (not included)
    
    %% 8. Synovial Reaction (not included)
    
    %% 9. Medial Joint Capsule Repair and 10. Growth Plate Thickness (both medial and lateral)
    
    SYN = sqrt((sx1-sx2)^2+(sy1-sy2)^2);
    
    medGP = sqrt((mgpx1-mgpx2)^2+(mgpy1-mgpy2)^2);
    latGP = sqrt((lgpx1-lgpx2)^2+(lgpy1-lgpy2)^2);
    
    
    tend = fix(clock);
    
    if tstart(4)<10
        hstart = ['0' int2str(tstart(4))];
    else
        hstart = [int2str(tstart(4))];
    end
    
    if tstart(5)<10
        mstart= ['0' int2str(tstart(5))];
    else
        mstart = [int2str(tstart(5))];
    end
    
    if tstart(6)<10
        sstart= ['0' int2str(tstart(6))];
    else
        sstart = [int2str(tstart(6))];
    end
    
    
    
    if tend(4)<10
        hend = ['0' int2str(tend(4))];
    else
        hend = [int2str(tend(4))];
    end
    
    if tend(5)<10
        mend = ['0' int2str(tend(5))];
    else
        mend = [int2str(tend(5))];
    end
    
    if tend(6)<10
        send = ['0' int2str(tend(6))];
    else
        send = [int2str(tend(6))];
    end
    
    
    
    
    
    
    
    masterhist{1,1} = 'Sample ID';
    masterhist{1,2} = 'Medial Compartment Tibial Plateau Width (pixels)';
    masterhist{1,3} = 'Cartilage Matrix Loss Width 0% Lesion Depth (pixels)';... cartilage matrix loss width 0% depth, pixels
        masterhist{1,4} = 'Cartilage Matrix Loss Width 50% Lesion Depth (pixels)';... cartilage matrix loss width 50% depth, pixels
        masterhist{1,5} = 'Cartilage Matrix Loss Width 1000% Lesion Depth (pixels)';... cartilage matrix loss width 100% depth, pixels
        masterhist{1,6} = 'Zone 1 (Medial) Cartilage Degneration Score';... cartilage degeneration score zone1
        masterhist{1,7} = 'Zone 2 (Central) Cartilage Degneration Score';... cartilage degeneration score zone2
        masterhist{1,8} = 'Zone 3 (Lateral) Cartilage Degneration Score';... cartilage degeneration score zone3
        masterhist{1,9} = 'Total Cartilage Degneration Width (% of Tibial Plateau)';... total cartilage degeneration width, percentage plateau width
        masterhist{1,10} = 'Significant Cartilage Degneration Width (% of Tibial Plateau)';... significant cartilage degeneration width, percentage plateau width
        masterhist{1,11} = 'Zone 1 (Medial) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 1
        masterhist{1,12} = 'Zone 2 (Central) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 2
        masterhist{1,13} = 'Zone 3 (Lateral) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 3
        masterhist{1,14} = 'Osteophyte Size (pixels)';... osteophyte length, pixels
        masterhist{1,15} = 'Medial Joint Capsule Repair (pixels)';... medial joint capsule repail, pixels
        masterhist{1,16} = 'Medial Compartment Medial Growth Plate Thickness (pixels)';... medial growth plate thickness, pixels
        masterhist{1,17} = 'Medial Compartment Lateral Growth Plate Thickness (pixels)';... lateral growth plate thickness, pixels
        masterhist{1,18} = '# pixel rows in image';
    masterhist{1,19} = '# pixel columns in image';
%     masterhist{1,20} = 'start timestamp';
%     masterhist{1,21} = 'end timestamp';
    
    
    
    masterhist{p+1,1} = ID; %'Sample ID';
    masterhist{p+1,2} = plateau; %'Medial Compartment Tibial Plateau Width (pixels)';
    masterhist{p+1,3} = width0; %'Cartilage Matrix Loss Width 0% Lesion Depth (pixels)';... cartilage matrix loss width 0% depth, pixels
    masterhist{p+1,4} = width50; %'Cartilage Matrix Loss Width 50% Lesion Depth (pixels)';... cartilage matrix loss width 50% depth, pixels
    masterhist{p+1,5} = width100; %'Cartilage Matrix Loss Width 1000% Lesion Depth (pixels)';... cartilage matrix loss width 100% depth, pixels
    masterhist{p+1,6} = CDSz1; %'Zone 1 (Medial) Cartilage Degneration Score';... cartilage degeneration score zone1
    masterhist{p+1,7} = CDSz2;%'Zone 2 (Central) Cartilage Degneration Score';... cartilage degeneration score zone2
    masterhist{p+1,8} = CDSz3;%'Zone 3 (Lateral) Cartilage Degneration Score';... cartilage degeneration score zone3
    masterhist{p+1,9} = 100*TCDW/plateau; %'Total Cartilage Degneration Width (% of Tibial Plateau)';... total cartilage degeneration width, percentage plateau width
    masterhist{p+1,10} = 100*SCDW/plateau; %'Significant Cartilage Degneration Width (% of Tibial Plateau)';... significant cartilage degeneration width, percentage plateau width
    masterhist{p+1,11} = 100*ZDRz1; %'Zone 1 (Medial) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 1
    masterhist{p+1,12} = 100*ZDRz2;%'Zone 2 (Central) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 2
    masterhist{p+1,13} = 100*ZDRz3;%'Zone 3 (Lateral) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 3
    masterhist{p+1,14} = OST; %'Osteophyte Size (pixels)';... osteophyte length, pixels
    masterhist{p+1,15} = SYN; %'Medial Joint Capsule Repair (pixels)';... medial joint capsule repail, pixels
    masterhist{p+1,16} = medGP; %'Medial Compartment Medial Growth Plate Thickness (pixels)';... medial growth plate thickness, pixels
    masterhist{p+1,17} = latGP;%'Medial Compartment Lateral Growth Plate Thickness (pixels)';... lateral growth plate thickness, pixels
    masterhist{p+1,18} = r; %'# pixel rows in image';
    masterhist{p+1,19} = c; %'# pixel columns in image';
%     masterhist{p+1,20} = [hstart ':' mstart ':' sstart]; %'start timestamp';
%     masterhist{p+1,21} = [hend ':' mend ':' send];%'end timestamp';
    
    
    
    t1 = subplot(3,4,5,'replace');
    text(0,0,sprintf('Complete.'), 'Parent', t1); axis off
    
    close all
    
    
    dr = 0;
    if exist(wfile)
        data = xlsread(wfile);
        [dr dc] = size(data);
        masterhist = masterhist(2:end,:);
        xlswrite(wfile,masterhist,1,['A' num2str(dr+2)])
    else
        xlswrite(wfile,masterhist,1,'A1')
    end
    
    
    
    
    
end
%
% clearvars( '-except', 'masterhist', 'p', 'pics', 'OARSI1', 'OARSI2', 'OARSI3', 'OARSI4', 'OARSI5', 'OARSI6', 'OARSI7', 'OARSI8', 'OARSI9', 'lines');
%
% dr = 0;
% if exist(wfile)
%     data = xlsread(wfile);
%     [dr dc] = size(data);
%     masterhist = masterhist(2:end,:);
%     xlswrite(wfile,masterhist,1,['A' num2str(dr+2)])
% else
%     xlswrite(wfile,masterhist,1,'A1')
% end
%


% xlswrite('GEKO_HEK.xls',masterhist,1,['A' num2str(dr+2)])
% xlswrite('GEKO_HEK.xls',masterhist(2,:),1,['A' num2str(OutputLine)])


display('GEKO has finished.')


% clear -EXCEPT masterhist









