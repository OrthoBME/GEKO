function GEKO_opensourcev6

% GEKO is an automated GUI to evaluate knee OA in frontally sectioned
% rodent histological images. It calculates 8 of the 10 OARSI recommended
% grades from the OARSI histopathology initiative (Gerwin 2010).

%% FUNCTION NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lines can be set to 'hide' or 'show' to have lines appear as you draw them
% or to suppress the lines

% wfile is the file name for summary file, should include extension. IE 'GEKO_Summary.xls'

% getptsNoDoubleClick was provided at 
% https://cbi.nyu.edu/svn/mrTools/trunk/mrUtilities/MatlabUtilities/
% by the NYU CBI (info@cbi.nyu.edu)

%% EDIT LOG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BYMJ 6/15/17: Added back button and disabled getpts double-click to end
% selection. Fixed so that images may be selected from any folder (does not
% have to be in same location as GEKO function). Fixed file name finding
% for general case. Added try loop and error catching. Cleaned up comments.
% Added comments. Adjusted calculations to account for getpts change.

% BYMJ 6/16/17: Added dialogue box questions to set initial parameters and 
% remove all print to command line information to prep for .exe packaging.
% Added zoom functionality.

%%
% Turns off warnings for the session - BYMJ
warning('off','all')
warning

% Creates question dialogue box to set lines mode. Hide setting is selected
% as default. - BYMJ
choice = questdlg('Show or hide output lines during grading?', ...
	'Set Lines', ...
	'show','hide','hide');
% Handle response
switch choice
    case 'show'
        lines = 'show';
    case 'hide'
        lines = 'hide';
end

% Requests user set a summary file name. -BYMJ
prompt = {'Enter desired summary file name'};
dlg_title = 'Summary title';
num_lines = 1;
defaultans = {'GEKO_Summary.xls'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
wfile = answer{1};

% Give users instruction to select image files. -BYMJ
waitfor(msgbox({'Please select images to grade.' 'Use shift + left mouse to select more than one file.',...
    'Files may only be selected from one directory at a time.'},'Select Image Help'));

% Select image files to work with
[fname1,pname1]=uigetfile({'*.*';'*.jpg';'*.png';'*.tiff'},'Select Image(s)','MultiSelect', 'on');

% This is needed for the multiselect.  Always put the file in this format.
% -BYMJ
fname1 = cellstr(fname1); 

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
Lablogo = imread('Lablogo.jpg');

% This for loop steps through each selected file. -BYMJ
for p = 1:length(fname1) % This for loop and everything before this line should be part of the GEKO shell code. -BYMJ
try % If a file crashes, continue on to the next file in "fname1". -BYMJ

    % Store current filename. -BYMJ
    filename = fullfile ( pname1, fname1{p} );
    % Break up filename into path, name, and extension.
    [~,name,~] = fileparts(filename);  
    
    % Read in the current file and set an identifier for the trial. -BYMJ
    I = imread(filename);
    ID = name;
    [r, c, ~] = size(I);   
    
    %% Create the marking points
    
    % Variable to indicate whether a regrade is needed. -BYMJ
    redo = 'y';
    while redo == 'y'      
        redo = 'n'; % Assume the image will not need to be regraded until the user inputs otherwise. -BYMJ
        
        % Create the Main GUI figure. -BYMJ
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on
        set(gcf,'Color','w');
        
        % Plot the reference image. -BYMJ
        subplot(3,4,1), image(OARSI1); axis off
        title('Reference Image')
        axis image 
        
        % Display the OrthoBME logo. -BYMJ
        subplot(3,4,9), image(Lablogo); axis off
        axis image
        subplot(3,4,[2:4,6:8,10:12]), image(I)
        ax = gca;
        title(ID,'interpreter','none')
        axis image
        
        linethick = 2;
           
        % Define tibial plateau
        clear px1 py1
               
        t1 = subplot(3,4,5, 'replace');
        fixfont = text(0,0,sprintf('Please click the \ntwo endpoints defining \nthe tibial plateau.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        [px1, py1] = getptsNoDoubleClick(ax);
        hold on       
        if exist('px1','var')
            subplot(3,4,[2:4,6:8,10:12]), plot(px1,py1,'ko');
        end
        
        % Store x and y coordinates in a structure to save later.
        rawpoints.tibplateauX = px1;
        rawpoints.tibplateauY = py1; 
        
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            subplot(3,4,[2:4,6:8,10:12]), plat = plot(px1,py1,'k');
            set(plat,'LineWidth',linethick);
        end
        
        % Define joint capsule thickness
        clear sx1 sy1
        subplot(3,4,1), image(OARSI2), axis off
        title('Reference Image')
        axis image         
        
        t1 = subplot(3,4,5,'replace');
        fixfont = text(0,0,sprintf('Please click the \ntwo endpoints defining \nthe joint capsule.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        [sx1, sy1] = getptsNoDoubleClick(ax);
        hold on
        if exist('sx1','var')
            subplot(3,4,[2:4,6:8,10:12]), plot(sx1,sy1,'go');
        end
        
        rawpoints.jointcapsuleX = sx1;
        rawpoints.jointcapsuleY = sy1;         
        
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            subplot(3,4,[2:4,6:8,10:12]), plat = plot(sx1,sy1,'g');
            set(plat,'LineWidth',linethick);
        end
        
        % Define growth plate medial and lateral thickness
        clear mgpx1 mgpy1 lgpx1 lgpy1
        subplot(3,4,1), image(OARSI3), axis off
        title('Reference Image')
        axis image  
        
        t1 = subplot(3,4,5,'replace');
        fixfont = text(0,0,sprintf('Please click the \ntwo endpoints defining the \nthickness of the medial half \nof the medial compartment \ngrowth plate.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        [mgpx1, mgpy1] = getptsNoDoubleClick(ax);
        hold on
        if exist('mgpx1','var')
            subplot(3,4,[2:4,6:8,10:12]), plot(mgpx1,mgpy1,'mo');
        end
        
        rawpoints.medgrowthplateX = mgpx1;
        rawpoints.medgrowthplateY = mgpy1; 
        
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            subplot(3,4,[2:4,6:8,10:12]), plat = plot(mgpx1,mgpy1,'m');
            set(plat,'LineWidth',linethick);
        end                
        
        subplot(3,4,1), image(OARSI4), axis off
        title('Reference Image')
        axis image    
        
        t1 = subplot(3,4,5,'replace');
        fixfont = text(0,0,sprintf('Please click the two \nendpoints defining the \nthickness of the lateral half \nof the medial compartment \ngrowth plate.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        [lgpx1, lgpy1] = getptsNoDoubleClick(ax);   
        hold on
        if exist('lgpx1','var')
            subplot(3,4,[2:4,6:8,10:12]), plot(lgpx1,lgpy1,'mo');
        end
        
        rawpoints.latgrowthplateX = lgpx1;
        rawpoints.latgrowthplateY = lgpy1;         
                
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            subplot(3,4,[2:4,6:8,10:12]), plat = plot(lgpx1,lgpy1,'m:');
            set(plat,'LineWidth',linethick);
        end        
        
        %Trace osteochondral interface
        clear ocx1 ocy1
        subplot(3,4,1), image(OARSI5), axis off
        title('Reference Image')
        axis image
        
        % Allows you to zoom in and out on the image. - BYMJ
        t1 = subplot(3,4,5, 'replace');
        fixfont = text(0,0,sprintf('Click and \ndrag an ROI \nto zoom in on \nif desired. \nEnter to continue.'),...
            'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        zoom on;
        pause() % you can zoom with your mouse and when your image is okay, you press any key. -BYMJ
        zoom off; % to escape the zoom mode -BYMJ     
        
        t1 = subplot(3,4,5,'replace');
        fixfont = text(0,0,sprintf('Please trace the \nosteochondral interface and \npress enter when done.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        clear ocy ocx1 ocy1

        [ocx1, ocy1] = getptsNoDoubleClick(ax);
        hold on
        if exist('ocx1','var')
            ocx = zeros(size(ocx1,1),1);
            ocy = zeros(size(ocy1,1),1);
        for i = 1:size(ocx1,1)
            ocx(i) = ocx1(i);
            ocy(i) = ocy1(i);
        end
            subplot(3,4,[2:4,6:8,10:12]), plot(ocx,ocy,'yo');
        end
        
        rawpoints.osteointerfaceX = ocx1;
        rawpoints.osteointerfaceY = ocy1; 
        
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            subplot(3,4,[2:4,6:8,10:12]), oc = plot(ocx,ocy,'y');
            set(oc,'LineWidth',linethick);
        end     
        zoom out;
        
        % Define total cartilage degeneration width
        clear tcdwx1 tcdwy1
        subplot(3,4,1), image(OARSI6), axis off
        title('Reference Image')
        axis image

        % Allows you to zoom in and out on the image. - BYMJ
        t1 = subplot(3,4,5, 'replace');
        fixfont = text(0,0,sprintf('Click and \ndrag an ROI \nto zoom in on \nif desired. \nEnter to continue.'),...
            'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        zoom on; % This will allow zoom once to degenerated ROI. Zoom out will occur after all lesion features are identified. -BYMJ
        pause() 
        zoom off; 
        
        t1 = subplot(3,4,5,'replace');
        fixfont = text(0,0,sprintf('If present, please click \nthe two endpoints defining \nthe total cartilage \ndegeneration width. \nIf not present, \njust press enter.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        [tcdwx1, tcdwy1] = getptsNoDoubleClick(ax);
        hold on
        if exist('tcdwx1','var')
            subplot(3,4,[2:4,6:8,10:12]), plot(tcdwx1,tcdwy1,'co');
        end
        
        rawpoints.totalcartdegenX = tcdwx1;
        rawpoints.totalcartdegenY = tcdwy1; 
        
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            if exist('tcdwx1','var')
                if size(tcdwx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot(tcdwx1,tcdwy1,'c');
                    set(plat,'LineWidth',linethick);
                end
            end
        end
        
        % Define significant cartilage degeneration width
        clear scdwx1 scdwy1 
        subplot(3,4,1), image(OARSI7), axis off
        title('Reference Image')
        axis image

        t1 = subplot(3,4,5,'replace');
        fixfont = text(0,0,sprintf('If present, please click \nthe two endpoints defining \nthe significant cartilage \ndegeneration width. \nIf not present, \njust press enter.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        [scdwx1, scdwy1] = getptsNoDoubleClick(ax);
        hold on
        if exist('scdwx1','var')
            subplot(3,4,[2:4,6:8,10:12]), plot(scdwx1,scdwy1,'bo');
        end
        
        rawpoints.sigcartdegenX = scdwx1;
        rawpoints.sigcartdegenY = scdwy1; 
        
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            if exist('scdwx1','var')
                if size(scdwx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot(scdwx1,scdwy1,'b:');
                    set(plat,'LineWidth',linethick);
                end
            end
        end
        
        % Define cartilage lesion
        subplot(3,4,1), image(OARSI8), axis off
        title('Reference Image')
        axis image
        clear cx cy cx1 cy1

        t1 = subplot(3,4,5,'replace');
        fixfont = text(0,0,sprintf('If present, trace \nthe cartilage lesion \nand press enter when done. \nIf not present, \njust press enter.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        [cx1, cy1] = getptsNoDoubleClick(ax);
        hold on
        if size(cx1,1)>0
            cx = zeros(size(cx1,1),1);
            cy = zeros(size(cy1,1),1);
            for i = 1:size(cx1,1)
                cx(i) = cx1(i);
                cy(i) = cy1(i);
            end
            subplot(3,4,[2:4,6:8,10:12]), plot(cx,cy,'ro');
        end
        
        rawpoints.lesionX = cx1;
        rawpoints.lesionY = cy1; 
        
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            if exist('cx','var')
                subplot(3,4,[2:4,6:8,10:12]), fill(cx,cy,'r');
            else
                subplot(3,4,[2:4,6:8,10:12]);
            end
        else
            subplot(3,4,[2:4,6:8,10:12]);
        end
        zoom out;
        
        %Define osteophyte thickness
        clear ostx1 osty1
        subplot(3,4,1), image(OARSI9), axis off
        title('Reference Image')
        axis image

        t1 = subplot(3,4,5,'replace');
        fixfont = text(0,0,sprintf('If present, please click \nthe two endpoints defining \nthe thickest diameter \nof the osteophyte \nfrom the chondro-osseous \njunction to the cartilage \nsurface. If not present, \njust press enter.'), 'Parent', t1); axis off
        fixfont.FontSize = 16;
        
        [ostx1, osty1] = getptsNoDoubleClick(ax);
        hold on
        if exist('ostx1','var')
            if size(ostx1)>0
                subplot(3,4,[2:4,6:8,10:12]), plot(ostx1,osty1,'yo');
            end
        end
        
        rawpoints.osteophyteX = ostx1;
        rawpoints.osteophyteY = osty1; 
        
        linecomp = strcmp(lines,'show');
        if linecomp == 1
            if exist('ostx1','var')
                subplot(3,4,[2:4,6:8,10:12]), plat = plot(ostx1,osty1,'y:');
                set(plat,'LineWidth',linethick);
            end
        end    
              
        % Plot all lines if the user wanted hidden lines during grading.
        linecomp = strcmp(lines,'show');
        if linecomp == 0
            hold on
            
            subplot(3,4,[2:4,6:8,10:12]), plat = plot(px1,py1,'k');
            set(plat,'LineWidth',linethick);            
            
            subplot(3,4,[2:4,6:8,10:12]), plat = plot(sx1,sy1,'g');
            set(plat,'LineWidth',linethick);            
            
            subplot(3,4,[2:4,6:8,10:12]), plat = plot(mgpx1,mgpy1,'m');
            set(plat,'LineWidth',linethick);
            
            subplot(3,4,[2:4,6:8,10:12]), plat = plot(lgpx1,lgpy1,'m:');
            set(plat,'LineWidth',linethick);            
            
            subplot(3,4,[2:4,6:8,10:12]), oc = plot(ocx,ocy,'y');
            set(oc,'LineWidth',linethick);            
            
            if exist('tcdwx1','var')
                if size(tcdwx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot(tcdwx1,tcdwy1,'c');
                    set(plat,'LineWidth',linethick);
                end
            end            
            
            if exist('scdwx1','var')
                if size(scdwx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot(scdwx1,scdwy1,'b:');
                    set(plat,'LineWidth',linethick);
                end
            end           
            
            if exist('cx','var')
                if size(cx)>0
                    subplot(3,4,[2:4,6:8,10:12]), fill(cx,cy,'r');
                end
            end            
            
            if exist('ostx1','var')
                if size(ostx1)>0
                    subplot(3,4,[2:4,6:8,10:12]), plat = plot(ostx1,osty1,'y:');
                    set(plat,'LineWidth',linethick);
                end
            end            
        end   
        
        % Creates question dialogue box to regrade image. No is selected
        % as default. - BYMJ
        choice = questdlg('Would you like to regrade the image?', ...
            'Regrade', ...
            'yes','no','no');
        % Handle response
        switch choice
            case 'yes'
                redo = 'y';
                clear rawpoints
            case 'no'
                redo = 'n';
        end          
    end  
%% Begin processing collected data

    % Create a waitbar so user knows the code is running. -BYMJ
    wait = waitbar(0,'Processing data. Please wait.','WindowStyle','modal');
    
    % Save graded image.
    saveas(gcf,fullfile(pname1, [ID '.png']), 'png'); 
    close
    
    waitbar(1/10); 
    
    % Save rawpoints structure and clear the variable. -BYMJ   
    save(ID,'rawpoints');
    clear rawpoints    
    
    waitbar(2/10);    
    
    % This would be the place to move calculations into a separate function
    % that loads the necessary variables and processes and saves the data.
    % -BYMJ
    
    %% 1. Cartilage Matrix Loss Width (0%, 50%, and 100% lesion depth)
    
    plateau = sqrt((px1(1)-px1(2))^2+(py1(1)-py1(2))^2);
    theta = atand((py1(1)-py1(2))/(px1(1)-px1(2)));
    
    if exist('cx','var')
        
        %Determining whether tibial line (black) is above or below lesion)
        m = (py1(2)-py1(1))/(px1(2)-px1(1));
        b = py1(1)-m*px1(1);
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
        
        h1 = sqrt((px1(1)-cx(1))^2+(py1(1)-cy(1))^2);
        h2 = sqrt((px1(2)-cx(1))^2+(py1(2)-cy(1))^2);
        h3 = sqrt((px1(1)-px1(2))^2+(py1(1)-py1(2))^2);
        
        s = (h1+h2+h3)/2;
        H = sqrt((4*s*(s-h1)*(s-h2)*(s-h3))/h3^2);
        
        dely = tibline*H*sind(90-theta);
        delx = tibline*H*cosd(90-theta);
        
        morey = tibline*round(r/30)*sind(90-theta);
        morex = tibline*round(r/30)*cosd(90-theta);
        
        
%       oldplat = [px1(1) px1(2) py1(1) py1(2)];
        px = px1;
        py = py1;
        px1 = (px(1)-delx)+(tibline)*morex;
        px2 = (px(2)-delx)+(tibline)*morex;
        py1 = (py(1)+dely)+(tibline*-1)*morey;
        py2 = (py(2)+dely)+(tibline*-1)*morey;       
                        
        % Creating the lesion on the image
        Ltemp = zeros(r,c);
        for i = 2:size(cx,1)
            rpts = linspace(cy(i-1),cy(i),1000);   %# A set of row points for the line
            cpts = linspace(cx(i-1),cx(i),1000);   %# A set of column points for the line
            index = sub2ind([r c],round(rpts),round(cpts));
            Ltemp(index) = 1;
            if i == size(cx,1)
                rpts = linspace(cy(1),cy(i),1000);   %# A set of row points for the line
                cpts = linspace(cx(1),cx(i),1000);   %# A set of column points for the line
                index = sub2ind([r c],round(rpts),round(cpts));
                Ltemp(index) = 1;
            end
        end
        
        for i = 2:size(ocx,1)
            rpts = linspace(ocy(i-1),ocy(i),1000);   %# A set of row points for the line
            cpts = linspace(ocx(i-1),ocx(i),1000);   %# A set of column points for the line
            index = sub2ind([r c],round(rpts),round(cpts));
            Ltemp(index) = 1;
        end
        
        stats = regionprops(Ltemp,'FilledImage','BoundingBox','Centroid','MajorAxisLength','Extrema', 'ConvexArea', 'ConvexHull','ConvexImage','Image','Orientation');
        lesion = imrotate(stats.FilledImage,theta);
        
        lesion = imdilate( lesion,strel('square',3));
        lesion =imerode(lesion,strel('square',2));
        [lesion, lnum] = bwlabel(lesion(:,:,1));
        
        
        if lnum>1
            index = [min(find(sum(lesion==1,1))) min(find(sum(lesion==1,2)))];
            target = find(index==max(index));
            lesion = lesion(:,min(find(sum(lesion==target,1))):max(find(sum(lesion==target,1))));
            
            [csr, csc] = find(lesion==target);
            cart_surf = min(csr);
            switch target
                case 1
                    [tmr, tmc] = find(lesion==2);
                    tidemark = min(tmr);
                case 2
                    [tmr, tmc] = find(lesion==1);
                    tidemark = min(tmr);
            end
            lesion = lesion(cart_surf:tidemark,:);
            lesion = lesion==target;
            
        else
            Ltemp = zeros(r,c);
            for i = 2:size(cx,1) % Changed to size(cx,1) from size(cx,2) - DFX
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
            [lesion, lnum] = bwlabel(lesion(:,:,1));
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
        
        % corrects the rotation angle if the lesion does not quite match
        % the tibial plateau angle
        alpha = atand((left(2)-right(2))/(left(1)-right(1)));
        lesion = imrotate(lesion,.8*alpha);
        lesion = lesion(min(find(sum(lesion,2)>0)):end,:); %cut off the blank rows at the top
        
        lesionbrackets = linspace(1,size(lesion,1),12);
        
        lesion0 = lesion(1:max([round(lesionbrackets(2)) round(left(2))+1]),:);
        lesion50 = lesion(round(lesionbrackets(6)):round(lesionbrackets(7)),:);
        lesion100 = lesion(round(lesionbrackets(11)):end,:);
                
        %two approaches via solid rows and absolute columns:
        %widest row of lesion
        width0 = max(sum(lesion0,2));
        width50 = max(sum(lesion50,2));
        width100 = max(sum(lesion100,2));
        
        waitbar(3/10);
        
        %% Cartilage Degeneration Score
        lesfill = stats.FilledImage;
        Ltemp = zeros(r,c);
        Ltemp(round(stats.BoundingBox(2)):round(stats.BoundingBox(2))+size(lesfill,1)-1,...
            round(stats.BoundingBox(1)):round(stats.BoundingBox(1))+size(lesfill,2)-1)=lesfill;
        
        rpts = linspace(py1,py2,1000);
        cpts = linspace(px1,px2,1000);   %# A set of column points for the line
        index = sub2ind([r c],round(rpts),round(cpts));
        Ltemp(index) = 1;
        
        L = Ltemp;
        L = imdilate( L,strel('square',10));
        L =imerode(L,strel('square',8));
        [L, lnum] = bwlabel(L(:,:,1));
        
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
        
        waitbar(4/10);
                                    
%% Zonal Depth Ratio
        
        Ltemp(round(stats.BoundingBox(2)):round(stats.BoundingBox(2))+size(lesfill,1)-1,round(stats.BoundingBox(1)):round(stats.BoundingBox(1))+size(lesfill,2)-1)=lesfill;
        
        rpts = linspace(py1,py2,1000);   %# A set of row points for the line
        cpts = linspace(px1,px2,1000);   %# A set of column points for the line
        index = sub2ind([r c],round(rpts),round(cpts));
        Ltemp(index) = 1;
        
        for i = 2:size(ocx,2)
            rpts = linspace(ocy(i-1),ocy(i),1000);   %# A set of row points for the line
            cpts = linspace(ocx(i-1),ocx(i),1000);   %# A set of column points for the line
            index = sub2ind([r c],round(rpts),round(cpts));
            Ltemp(index) = 1;
        end
        
        L = Ltemp;
        L = imdilate( L,strel('square',8));
        L =imerode(L,strel('square',6));
        [L, ~] = bwlabel(L(:,:,1));        
        
        L = imrotate(L,theta);
        L = L(min(find(sum(L,2)>0)):max(find(sum(L,2)>0)),min(find(sum(L,1)>0)):max(find(sum(L,1)>0)));
        L = L(6:end,:);
        L = L(min(find(sum(L,2)>0)):max(find(sum(L,2)>0)),min(find(sum(L,1)>0)):max(find(sum(L,1)>0)));
                
        [L, lnum] = bwlabel(L);
        zones = linspace(1,size(L,2),7);
        
        zone1 = L(:,round(zones(5)):end);
        zone2 = L(:,round(zones(3)):round(zones(5)));
        zone3 = L(:,1:round(zones(3)));
        
        center = [zone1(:,round(zones(2))) zone2(:,round(zones(2))) zone3(:,round(zones(2)))];
        
        if lnum==1
            for i = 1:3
                if sum(center(:,i)==1)>0  % red=2, yellow=1
                    [holdmat, hnum] = bwlabel(center(:,i));
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
    
    waitbar(5/10); 
    
%%   Total and  Significant and total Cartilage Degeneration Width
    
    if exist('tcdwx1','var')
        if size(tcdwx1)>0
            TCDW = sqrt((tcdwx1(1)-tcdwx1(2))^2+(tcdwy1(1)-tcdwy1(2))^2);
        else
            TCDW = 0;
        end
    else
        TCDW = 0;
    end
    if exist('scdwx1','var')
        if size(scdwx1)>0
            SCDW = sqrt((scdwx1(1)-scdwx1(2))^2+(scdwy1(1)-scdwy1(2))^2);
        else
            SCDW = 0;
        end
    else
        SCDW = 0;
    end
    
    waitbar(6/10);
    
%% 6. Osteophyte    
    
    if exist('ostx1','var')
        if size(ostx1)>0
            OST = sqrt((ostx1(1)-ostx1(2))^2+(osty1(1)-osty1(2))^2);
        else
            OST = 0;
        end
    else
        OST = 0;
    end
    
    waitbar(7/10);    
    
    %% 7. Calcified Cartilage and Subchondral bone Damage Score (not included)
    
    %% 8. Synovial Reaction (not included)
    
    %% 9. Medial Joint Capsule Repair and 10. Growth Plate Thickness (both medial and lateral)
    
    SYN = sqrt((sx1(1)-sx1(2))^2+(sy1(1)-sy1(2))^2);
    
    medGP = sqrt((mgpx1(1)-mgpx1(2))^2+(mgpy1(1)-mgpy1(2))^2);
    latGP = sqrt((lgpx1(1)-lgpx1(2))^2+(lgpy1(1)-lgpy1(2))^2);
    
    waitbar(8/10);
        
%% Create master file labels    
    masterhist{1,1} = 'Sample ID';
    masterhist{1,2} = 'Medial Compartment Tibial Plateau Width (pixels)';
    masterhist{1,3} = 'Cartilage Matrix Loss Width 0% Lesion Depth (pixels)';... cartilage matrix loss width 0% depth, pixels
    masterhist{1,4} = 'Cartilage Matrix Loss Width 50% Lesion Depth (pixels)';... cartilage matrix loss width 50% depth, pixels
    masterhist{1,5} = 'Cartilage Matrix Loss Width 1000% Lesion Depth (pixels)';... cartilage matrix loss width 100% depth, pixels

% These are ordinal scores GEKO calculates, but have not been included in
% the final output. Can be included if desired.
%     masterhist{1,6} = 'Zone 1 (Medial) Cartilage Degneration Score';... cartilage degeneration score zone1
%     masterhist{1,7} = 'Zone 2 (Central) Cartilage Degneration Score';... cartilage degeneration score zone2
%     masterhist{1,8} = 'Zone 3 (Lateral) Cartilage Degneration Score';... cartilage degeneration score zone3

    masterhist{1,6} = 'Total Cartilage Degneration Width (% of Tibial Plateau)';... total cartilage degeneration width, percentage plateau width
    masterhist{1,7} = 'Significant Cartilage Degneration Width (% of Tibial Plateau)';... significant cartilage degeneration width, percentage plateau width
    masterhist{1,8} = 'Zone 1 (Medial) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 1
    masterhist{1,9} = 'Zone 2 (Central) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 2
    masterhist{1,10} = 'Zone 3 (Lateral) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 3
    masterhist{1,11} = 'Osteophyte Size (pixels)';... osteophyte length, pixels
    masterhist{1,12} = 'Medial Joint Capsule Repair (pixels)';... medial joint capsule repair, pixels
    masterhist{1,13} = 'Medial Compartment Medial Growth Plate Thickness (pixels)';... medial growth plate thickness, pixels
    masterhist{1,14} = 'Medial Compartment Lateral Growth Plate Thickness (pixels)';... lateral growth plate thickness, pixels
    masterhist{1,15} = '# pixel rows in image';
    masterhist{1,16} = '# pixel columns in image';

    waitbar(9/10);
    
    masterhist{p+1,1} = ID; %'Sample ID';
    masterhist{p+1,2} = plateau; %'Medial Compartment Tibial Plateau Width (pixels)';
    masterhist{p+1,3} = width0; %'Cartilage Matrix Loss Width 0% Lesion Depth (pixels)';... cartilage matrix loss width 0% depth, pixels
    masterhist{p+1,4} = width50; %'Cartilage Matrix Loss Width 50% Lesion Depth (pixels)';... cartilage matrix loss width 50% depth, pixels
    masterhist{p+1,5} = width100; %'Cartilage Matrix Loss Width 100% Lesion Depth (pixels)';... cartilage matrix loss width 100% depth, pixels
    
    masterhist{p+1,6} = 100*TCDW/plateau; %'Total Cartilage Degneration Width (% of Tibial Plateau)';... total cartilage degeneration width, percentage plateau width
    masterhist{p+1,7} = 100*SCDW/plateau; %'Significant Cartilage Degneration Width (% of Tibial Plateau)';... significant cartilage degeneration width, percentage plateau width
    masterhist{p+1,8} = 100*ZDRz1; %'Zone 1 (Medial) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 1
    masterhist{p+1,9} = 100*ZDRz2;%'Zone 2 (Central) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 2
    masterhist{p+1,10} = 100*ZDRz3;%'Zone 3 (Lateral) Zonal Depth Ratio of Lesion (% Full Cartilage Thickness)';... zonal depth ratio of lesions, zone 3
    masterhist{p+1,11} = OST; %'Osteophyte Size (pixels)';... osteophyte length, pixels
    masterhist{p+1,12} = SYN; %'Medial Joint Capsule Repair (pixels)';... medial joint capsule repail, pixels
    masterhist{p+1,13} = medGP; %'Medial Compartment Medial Growth Plate Thickness (pixels)';... medial growth plate thickness, pixels
    masterhist{p+1,14} = latGP;%'Medial Compartment Lateral Growth Plate Thickness (pixels)';... lateral growth plate thickness, pixels
    masterhist{p+1,15} = r; %'# pixel rows in image';
    masterhist{p+1,16} = c; %'# pixel columns in image';
    
        
    waitbar(10/10);
    close(wait);    
       
if length(fname1) == 1
    break
end    
catch MException
    if exist('wait','var') == 1
        close(wait);
    end
    waitfor(msgbox('ERROR encountered, continue to next image.','ERROR','WindowStyle','modal'))
    beep
end    
end

% Does not try to save if the trial crashed.
if exist('masterhist','var') == 1
    % Creates a dialogue box indicated excel sheet is saving. Indicates when
    % code is complete. -BYMJ
    d = msgbox('Saving data to excel file...','Save','WindowStyle','modal');
    set(findobj(d,'style','pushbutton'),'Visible','off')
    
    system('taskkill /F /IM EXCEL.EXE');
    if exist(fullfile(pname1,wfile),'file')
        data = xlsread(fullfile(pname1,wfile));
        system('taskkill /F /IM EXCEL.EXE');
        [dr, ~] = size(data);
        masterhist = masterhist(2:end,:);
        xlswrite(fullfile(pname1,wfile),masterhist,1,['A' num2str(dr+2)])
    else
        xlswrite(fullfile(pname1,wfile),masterhist,1,'A1')
    end

    set(findobj(d,'Tag','MessageBox'),'String','GEKO has finished.')
    set(findobj(d,'style','pushbutton'),'Visible','on')
    close all
end
end












