stims = {'afasa1.jpg', 'afasa2.jpg', 'afasa3.jpg', 'afasa4.jpg'};

for im = 1:length(stims)
    x=imread(stims{im});
    y = 255 - x;
    imwrite(y,[stims{im}(1:6), '_inv', stims{im}(7:end)]);
end