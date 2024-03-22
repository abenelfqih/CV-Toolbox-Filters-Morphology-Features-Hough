function fftshow(f, type)

if nargin < 2
    type = 'log';
end

if (strcmp(type, 'log'))
    fl = log(1 + abs(f));
    fm = max(fl(:));
    imshow(im2uint8(fl / fm))
elseif (strcmp(type, 'abs'))
    fa = abs(f);
    fm = max(fa(:));
    imshow(fa / fm)
else
    error('TYPE doit ?tre de type "abs" ou "log"');
end