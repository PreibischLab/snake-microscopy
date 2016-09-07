img = imread('a.png');
img = rgb2gray(img);

figure; imshow(img); hold on;
for i=1:14
    load([num2str(i) '.mat']);
    plot([Xs;Xs], [Ys;Ys], 'LineWidth', 2)
end

Xs = nan;
Ys = nan;

[Xs, Ys] = setManualSnake(img, 12);
 [Xs, Ys] = resampleSnake(img, Xs, Ys);
save('14','Xs','Ys');