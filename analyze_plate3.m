function analyze_plate3
%%%
%   V.3
%   incl. gaussian background correction
%   for fixed and stained plates

%---------
close all;
clear all;
format compact;
%---------

%read img
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile('*.tif', 'choose plate image file ...');

mkdir(PATHNAME,FILENAME(1:length(FILENAME)-4));
NEWPATH=[PATHNAME,FILENAME(1:length(FILENAME)-4),'\'];

plate=imread(strcat(PATHNAME,FILENAME));
[rows,cols]=size(plate);

%top plate
platetop=plate;
platetop=rot90(platetop,3);
platetop02=imresize(platetop,0.2);
platetop=imcomplement(platetop);
platetop02bin=platetop02>35000;

%display top plate
figure;
subplot(2,1,1),
imshow(platetop02,[]);
subplot(2,1,2), imshow(platetop02bin,[]);

%find top plate wells
[centers,radii] = imfindcircles(platetop02bin,[160 200],'Sensitivity',...
    0.9500,'EdgeThreshold',0.30,'Method','PhaseCode',...
    'ObjectPolarity','Bright');
viscircles(centers,radii,'EdgeColor','r');

%combine data
position=cat(2,centers.*5,radii.*5);
%sort wells in Y direction
[values, order] = sort(position(:,2));
sortpos = position(order,:);
%sort wells in X direction
row=mat2cell(sortpos,[3 3],[3]);
row1=row{1,1};
row2=row{2,1};
[values, order] = sort(row1(1:3,1));
row1sort = row1(order,1:3);
[values, order] = sort(row2(1:3,1));
row2sort = row2(order,1:3);
sorted=cat(1,row1sort,row2sort);
%crop wells
figure;
array_to_export=[];
for i=1:6
    xmin=min(max((sorted(i,1)-floor((sorted(i,3)+5))),0),cols-1);
    ymin=min(max((sorted(i,2)-floor((sorted(i,3)+5))),0),rows-1);
    width=(floor(sorted(i,3))*2+10);
    height=(floor(sorted(i,3))*2+10);
    wellgray=imcrop(platetop,[xmin ymin width-1 height-1]);
    
    %mask wells
    [x,y]=meshgrid(-floor(-width/2:width/2-1),-floor(-width/2:width/2-1));
    c_mask=uint16((x.^2+y.^2)<=(sorted(i,3)-150)^2);
    mask=size(c_mask);  
    org_wellgray=wellgray.*c_mask;
    bwgray=wellgray;
    wellgray=wellgray.*c_mask;
    
    well_radius(i)=floor(sorted(i,3))-150;
    

%   save well image
    filename=['well_',num2str(i),'.tif'];
    imwrite(uint16(org_wellgray),[NEWPATH,filename],'TIFF');
    
 %background correction
    bckgrnd=imgaussfilt(bwgray,100);
   
    wellgray=wellgray-bckgrnd;

    %   save well image
    filename=['bkgrnd_',num2str(i),'.tif'];
    imwrite(uint16(wellgray),[NEWPATH,filename],'TIFF');    
  %   save well image
    filename=['bkgrnd2_',num2str(i),'.tif'];
    imwrite(uint16(bckgrnd),[NEWPATH,filename],'TIFF');    
      
%    
 
I_dilate = segment_well(uint16(wellgray),0.12);
I_dilate=I_dilate>0;
%figure;imshow(I_dilate,[]);

%   save dil image
    filename=['dil_',num2str(i),'.tif'];
    imwrite(uint16(I_dilate),[NEWPATH,filename],'TIFF');
    

    [o_mask,results]=regiongrowing(I_dilate,uint16(wellgray),0.1,75,0.120);
    
wellgray=double(wellgray);
overlay=wellgray.*o_mask;
filename=['overlay_',num2str(i),'.tif'];
imwrite(uint16(overlay),[NEWPATH,filename],'TIFF');

org_wellgray=im2double(org_wellgray);
overlay2=cat(3,org_wellgray,org_wellgray.*(1-overlay),org_wellgray.*(1-overlay));

filename=['overlay2_',num2str(i),'.tif'];
imwrite((overlay2),[NEWPATH,filename],'TIFF');

    for j=1:size(results,1)
        array_to_export = cat(1,array_to_export,[i results(j,1) results(j,2) results(j,3) results(j,4) results(j,5)]);
    end
    
    %plot wells
    hold on
    subplot(2,3,i);imshow(I_dilate,[]);
    title(i);
    hold off
end

array_to_export=cat(1,[well_radius(1) well_radius(2) well_radius(3) well_radius(4) well_radius(5) well_radius(6)], array_to_export);
tablename=[FILENAME(1:length(FILENAME)-4),'_results.csv'];
    csvwrite([NEWPATH,tablename],array_to_export);
end