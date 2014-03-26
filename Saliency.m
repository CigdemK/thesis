% clc, clear;
% img = initializeImage('../holywood2/getoutofcar/frames/video00148.jpg');
% params = defaultSaliencyParams;
% param.weights = [1 1 10]
% salmap = makeSaliencyMap(img,params);
% bigMap = salmap;
% bigMap.data = imresize(salmap.data,img.size(1:2));
% figure;
% displayMaps(salmap,1);

img = imread('../holywood2/getoutofcar/frames_original/video00148.jpg');
map = gbvs(img);
map_itti = ittikochmap(img);
figure;
show_imgnmap( img , map );
figure;
show_imgnmap( img , map_itti );