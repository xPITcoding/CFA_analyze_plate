function [Idilate]=segment_well(I,threshold)

%---------
%close all;
%clear all;
%---------
I=im2double(I);
%I=imread('well_1.tif');
Ibin=imbinarize(I,threshold);
closefac=strel('square',7);
Idilate=imdilate(Ibin,closefac);
openfac=strel('disk',3);
Idilate=imerode(Idilate,openfac);

%subplot(1,2,1), imshow(Ibin);
%subplot(1,2,2), imshow(Idilate);

end