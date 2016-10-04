img = imread('building_test_1.jpg');
img_copy = img;

figure,imshow(img_copy);

%img = uint8(bilateral_filter(double(img),5,5,15));
%figure,imshow(edge(rgb2gray(img_copy),'canny'));


% Do something here to make the windows stand out. Some preprocessing
% required. Increase contrast in image?

figure,imshow(img);

%F = fft2(img);

%F = fftshift(F); % Center FFT

%store_fft = F;
%figure,imshow(mat2gray(log(abs(store_fft)+1)));

%F = abs(F); % Get the magnitude
%F = log(F+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
%F = mat2gray(F); % Use mat2gray to scale the image between 0 and 1

%f1 = F(:,:,1);
%f2 = F(:,:,2);
%f3 = F(:,:,3);

%lim1 = (mean(f1(:)) + 0.8*max(f1(:)))/1.8;
%lim2 = (mean(f2(:)) + 0.8*max(f2(:)))/1.8;
%lim3 = (mean(f3(:)) + 0.8*max(f3(:)))/1.8;

%for i=1:size(F,1)
%    for j=1:size(F,2)
%        if F(i,j,1) < lim1
%            store_fft(i,j,1) = 0.0;
%        end
%        if F(i,j,2) < lim2
%            store_fft(i,j,2) = 0.0;
%        end
%        if F(i,j,3) < lim3
%            store_fft(i,j,3) = 0.0;
%        end
%    end
%end

%figure,imshow(mat2gray(log(abs(store_fft)+1)));

%final_img = uint8(ifft2(ifftshift(store_fft)));
%final_img = imgaussfilt(final_img,2);

%figure,imshow(final_img);

[level, EM] = graythresh(img);
BW = im2bw(img,level);

BW = imcomplement(BW);

BW = imerode(BW,strel('rectangle',[10 10]));
BW = imdilate(BW,strel('rectangle',[10 10]));

figure,imshow(BW);

conn = conndef(ndims(BW), 'minimal');
CC = bwconncomp(BW, conn);
D1 = regionprops(CC, 'area', 'perimeter'); 

areas = [D1.Area];
[meanArea] = mean(areas);
disp(meanArea);

L = labelmatrix(CC);

bwfinal = ismember(L, find((abs([D1.Area] - ([D1.Perimeter] ./ 4) .* ([D1.Perimeter] ./ 4)) <= [D1.Area] ./ 3) & ([D1.Area] >= meanArea/5.0)));

figure, imshow(bwfinal);

edges = imdilate(bwfinal,strel('rectangle',[5 5])) - imerode(bwfinal,strel('rectangle',[5 5]));

figure,imshow(edges);

total_img = img_copy;

for i=1:size(total_img,1)
    for j=1:size(total_img,2)
        if edges(i,j) ~= 0.0
            total_img(i,j,1) = 255.0;
            total_img(i,j,2) = 0.0;
            total_img(i,j,3) = 0.0;
        end
    end
end

figure, imshow(total_img);
