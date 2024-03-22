function varargout = menu1(varargin)
% MENU1 MATLAB code for menu1.fig
%      MENU1, by itself, creates a new MENU1 or raises the existing
%      singleton*.
%
%      H = MENU1 returns the handle to a new MENU1 or the handle to
%      the existing singleton*.
%
%      MENU1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MENU1.M with the given input arguments.
%
%      MENU1('Property','Value',...) creates a new MENU1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before menu1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to menu1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help menu1

% Last Modified by GUIDE v2.5 27-Feb-2024 21:15:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @menu1_OpeningFcn, ...
                   'gui_OutputFcn',  @menu1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before menu1 is made visible.
function menu1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to menu1 (see VARARGIN)

% Choose default command line output for menu1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes menu1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = menu1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function MenuFichier_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFichier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuBruit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBruit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuGaussien_Callback(hObject, eventdata, handles)
% hObject    handle to MenuGaussien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
% V?rifier si axes2 contient quelque chose
 
Imgb =imnoise(imageO,'gaussian', 0,0.01);
%Imgb=handles.Imgb;
axes(handles.axes2)
imshow(Imgb);
handles.courant_data=Imgb;
%subImage(handles.Imgb);
handles.output=hObject;

% --------------------------------------------------------------------
function menuPVetSel_Callback(hObject, eventdata, handles)
% hObject    handle to menuPVetSel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
% V?rifier si axes2 contient quelque chose

Imgb =imnoise(imageO,'Salt & Pepper', 0.01);
axes(handles.axes2)
imshow(Imgb);
handles.courant_data=Imgb;
%subImage(handles.Imgb);
handles.output=hObject;

% --------------------------------------------------------------------
function MenuOuvrir_Callback(hObject, eventdata, handles)
% hObject    handle to MenuOuvrir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.*');
handles.ima = imread(sprintf('%s',path,file));
axes(handles.axes1)
handles.courant_data = handles.ima;
subimage(handles.courant_data);
axes(handles.axes2)
% handles.ima_traite = 0;
% subimage(handles.ima_traite);

subimage(handles.courant_data);

handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function MenuEnregistrer_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEnregistrer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[file,path] = uiputfile('*.png','Enregistrer Votre Image ...');
imwrite(image, sprintf('%s',path,file),'png');

% --------------------------------------------------------------------
function MenuQuitter_Callback(hObject, eventdata, handles)
% hObject    handle to MenuQuitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --------------------------------------------------------------------
function menu_Transformation_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Transformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_filtrePassBas_Callback(hObject, eventdata, handles)
% hObject    handle to menu_filtrePassBas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_filtrePassHaut_Callback(hObject, eventdata, handles)
% hObject    handle to menu_filtrePassHaut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_filtreGradient_Callback(hObject, eventdata, handles)
% hObject    handle to menu_filtreGradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
image = double(image);
[m,n] = size(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
maskhor = [0,0,0;-1,0,1;0,0,0]; 
maskver = [0,-1,0;0,0,0;0,1,0];
for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3         
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);            
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
 end
%mymin=min(min(output))
%mymax=max(max(output))
for i=4:(m-3)
    for j=4:(n-3)       
        output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j));
    end 
end 
%outputhor=uint8(outputhor); 
%outputver=uint8(outputver); 
output=uint8(output); 

%b=uint8(b);
axes(handles.axes2);
imshow(output);

% --------------------------------------------------------------------
function menu_filtreRobert_Callback(hObject, eventdata, handles)
% hObject    handle to menu_filtreRobert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[n,m]=size(image);
image = double(image);
 %num = get(handles.slider1, 'value');
% set(handles.edit1, 'String', num2str(num));
for x=1:n-1
 for y=1:m-1
  b(x,y)= abs(uint8( double(image(x,y))-double(image(x+1,y+1))))+ abs(uint8( double(image(x,y+1)) - double(image(x+1,y))));
 end
end
    % num = get(handles.slider1, 'Value');
    % set(handles.txt1, 'String', num2str(num));
        %Seuillage
        [n,m]=size(image);
        for i=1:n-1
         for j=1:m-1
          if b(i,j) < 25
            b(i,j)=0;
          end
         end
        end
           %
  handles.ima_traite = b;
  axes(handles.axes2);
  imshow(b);
%Grrr
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_filtreLaplacienne_Callback(hObject, eventdata, handles)
% hObject    handle to menu_filtreLaplacienne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=handles.courant_data;
[n,m]=size(image);
image = double(image);
[n m]=size(image);
b=zeros(n,m);
M1=[-1 -1 -1;-1 8 -1;-1 -1 -1];
for i=2:n-1
    for j=2:m-1
        V=image((i-1:i+1),(j-1:j+1));
        S=V.*M1;
        b(i,j)=sum(S(:));
    end
end
b=uint8(b);
axes(handles.axes2);
     subimage(b);
% --------------------------------------------------------------------
function menu_filtre_gaussien_Callback(hObject, eventdata, handles)
% hObject    handle to menu_filtre_gaussien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/16)*[1 2 1 ;2 4 2 ; 1 2 1];
for x = 2 : n-1
    for y = 2 : m-1
    f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.axes2);
imshow(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_filtre_moyenneur_Callback(hObject, eventdata, handles)
% hObject    handle to menu_filtre_moyenneur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/9)*[1 1 1 ; 1 1 1 ; 1 1 1 ];
for x = 2 : n-1
    for y = 2 : m-1
      f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.axes2);
imshow(b);

% --------------------------------------------------------------------
function menu_filtre_mediane_Callback(hObject, eventdata, handles)
% hObject    handle to menu_filtre_mediane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
image=double(image);
[n,m]=size(image);
img=image;
for i=2:n-1
    for j=2:m-1
       fenetre=image(i-1:i+1,j-1:j+1);
       v=[fenetre(1,:) fenetre(2,:) fenetre(3,:)];
       sort(v);
       a=median(v);
       img(i,j)=a;
    end
end
b=uint8(img);
handles.ima_traite = b;
axes(handles.axes2);
imshow(b);
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Histogramme_Callback(hObject, eventdata, handles)
    % hObject    handle to Histogramme (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % R?cup?rer l'image depuis axes1
image = getimage(handles.axes1);

    

        d = length(size(image));
        if d==3
            I = rgb2gray(image);
        elseif d==2
            I = img
        end
axes(handles.axes1);
subimage(I);
[nl nc]=size(I);
v=double(I);
vec=[1:256];
l=0;
for k=0:255 
    for i=1:nl
        for j=1:nc
            if v(i,j)==k 
               l=l+1;
            end
        end
    end
    vec(k+1)=l;
    l=0;
end
axes(handles.axes2);plot(vec);





% --------------------------------------------------------------------
function menu_inverser_Callback(hObject, eventdata, handles)
% hObject    handle to menu_inverser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[l,c] = size(image);
v = image;
for i=1:l
    for j=1:c
        v(i,j) = -double(image(i,j))+255;
    end
end
v = uint8(v);
axes(handles.axes2);
imshow(v);


% --------------------------------------------------------------------
function menu_contrast_Callback(hObject, eventdata, handles)
% hObject    handle to menu_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
image = double(image);
[l,c] = size(image);
v = image;
for i=1:l
    for j=1:c
        fpixel = (image(i,j)-128)*5 +128;
        if (fpixel>255)
            fpixel=255;
        else
            if (fpixel<0)
                fpixel = 0;
            end
        end
        v(i,j) = fpixel;
    end
end
v = uint8(v);
axes(handles.axes2);
imshow(v);


% --------------------------------------------------------------------
function menu_liminosite_Callback(hObject, eventdata, handles)
% hObject    handle to menu_liminosite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[l,c] = size(image);
v = image;
for i=1:l
    for j=1:c
        pix = image(i,j)+50;
        if (pix>255)
            pix=255;
        else
            if(pix<0)
                pix = 0;
            end
        end
        v(i,j) = pix;
    end
end
v = uint8(v);
axes(handles.axes2);
imshow(v);


% --------------------------------------------------------------------
function Filtre_frequentiel_Callback(hObject, eventdata, handles)
% hObject    handle to Filtre_frequentiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function morphologie_Callback(hObject, eventdata, handles)
% hObject    handle to morphologie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PI_Callback(hObject, eventdata, handles)
% hObject    handle to PI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Susan_Callback(hObject, eventdata, handles)
% hObject    handle to Susan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Charger une image
im = handles.courant_data;

% V?rifier la dimension de l'image et la convertir en niveaux de gris si elle est en couleur
if length(size(im)) == 3
    image = double(rgb2gray(im));
elseif length(size(im)) == 2
    image = double(im);
else
    error('L''image doit ?tre en niveaux de gris ou en couleur.');
end

% =============================donn?es=====================================
rayon = 2;
alpha = 50;
r = 2;
alpha = alpha / 100;
[n, m] = size(image);

% ========================g?nerateur de mask===============================
mask = zeros(2*rayon + 1);
b = ones(rayon + 1);
for i = 1:rayon + 1
    for j = 1:rayon + 1
        if (rayon == 1)
            if (j > i)
                b(i, j) = 0;
            end
        else
            if (j > i + 1)
                b(i, j) = 0;
            end
        end
    end
end
mask(1:rayon + 1, rayon + 1:2*rayon + 1) = b;
mask(1:rayon + 1, 1:rayon + 1) = rot90(b);
mask0 = mask;
mask0 = flipdim(mask0, 1);
mask = mask0 + mask;
mask(rayon + 1, :) = mask(rayon + 1, :) - 1;

% ==========================r?ponse maximale===============================
max_reponse = sum(mask(:));

% =====================balayage de toute l'image===========================
f = zeros(n, m);
for i = (rayon + 1):n - rayon
    for j = (rayon + 1):m - rayon
        image_courant = image(i - rayon:i + rayon, j - rayon:j + rayon);
        image_courant_mask = image_courant .* mask;
        inteniste_cental = image_courant_mask(rayon + 1, rayon + 1);
        s = exp(-1 * (((image_courant_mask - inteniste_cental) / max_reponse).^6));
        somme = sum(s(:));
        %   si le centre du mask est un 0 il faut soustraire les zeros des filtres
        if (inteniste_cental == 0)
            somme = somme - sum(mask(:) == 0);
        end
        f(i, j) = somme;
    end
end

% =============selection et seuillage des points d'inter?t=================
ff = f(rayon + 1:n - (rayon + 1), rayon + 1:m - (rayon + 1));
minf = min(ff(:));
maxf = max(f(:));
fff = f;

d = 2*r + 1;
temp1 = round(n/d);
if (temp1 - n/d) < 0.5 && (temp1 - n/d) > 0
    temp1 = temp1 - 1;
end

temp2 = round(m/d);
if (temp2 - m/d) < 0.5 && (temp2 - m/d) > 0
    temp2 = temp2 - 1;
end

fff(n:temp1*d + d, m:temp2*d + d) = 0;

for i = (r + 1):d:temp1*d + d
    for j = (r + 1):d:temp2*d + d
        window = fff(i - r:i + r, j - r:j + r);
        window0 = window;
        [xx, yy] = find(window0 == 0);
        for k = 1:length(xx)
            window0(xx(k), yy(k)) = max(window0(:));
        end
        minwindow = min(window0(:));
        [y, x] = find(minwindow ~= window & window <= minf + alpha * (maxf - minf) & window > 0);
        [u, v] = find(minwindow == window);
        if length(u) > 1
            for l = 2:length(u)
                fff(i - r - 1 + u(l), j - r - 1 + v(l)) = 0;
            end
        end
        if ~isempty(x)
            for l = 1:length(y)
                fff(i - r - 1 + y(l), j - r - 1 + x(l)) = 0;
            end
        end
    end
end

seuil = minf + alpha * (maxf - minf);
[u, v] = find(minf <= fff & fff <= seuil);

% ==============affichage des resultats====================================
axes(handles.axes2);
imshow(im);
hold on;
plot(v, u, '.r', 'MarkerSize', 10);
nombre_de_point_dinteret = length(v);

% --------------------------------------------------------------------
function Harris_Callback(hObject, eventdata, handles)
% hObject    handle to Harris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.courant_data;
%==========================================================================
if(size(img,3)==3)
    display('l''image est en couleur')
    img=rgb2gray(img);
end
%==========================================================================

lambda=0.04;
sigma=1; seuil=200; r=6; w=5*sigma;
[m,n]=size(img);
imd=double(img);
dx=[-1 0 1
    -2 0 2
    -1 0 1]; % deriv?e horizontale : filtre de Sobel
dy=dx'; % deriv?e verticale : filtre de Sobel
g = fspecial('gaussian',max(1,fix(w)), sigma);
Ix=conv2(imd,dx,'same');
Iy=conv2(imd,dy,'same');
Ix2=conv2(Ix.^2, g, 'same');
Iy2=conv2(Iy.^2, g, 'same');
Ixy=conv2(Ix.*Iy, g,'same');
detM=Ix2.*Iy2-Ixy.^2;
trM=Ix2+Iy2;
R=detM-lambda*trM.^2;
%==========================================================================
R1=(1000/(1+max(max(R))))*R;
%==========================================================================          
[u,v]=find(R1<=seuil);
nb=length(u);
for k=1:nb
    R1(u(k),v(k))=0;
end
R11=zeros(m+2*r,n+2*r);
R11(r+1:m+r,r+1:n+r)=R1;
[m1,n1]=size(R11);
for i=r+1:m1-r
    for j=r+1:n1-r
        fenetre=R11(i-r:i+r,j-r:j+r);
        ma=max(max(fenetre));
        if fenetre(r+1,r+1)<ma
            R11(i,j)=0;
        end
    end
end
axes(handles.axes2);
imshow(img);
hold on;
R11=R11(r+1:m+r,r+1:n+r);
[x,y]=find(R11);
nb=length(x);
plot(y,x,'.r');
title('Detection des points d''interet');

% --------------------------------------------------------------------
function HarrisElectro_Callback(hObject, eventdata, handles)
% hObject    handle to HarrisElectro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img = handles.courant_data;

if size(img, 3) == 3
    img = rgb2gray(img);
end

lambda = 0.04;
sigma = 1; 
seuil = 100; 
r = 6; 
w = 5 * sigma;

[m, n] = size(img); 
imd = double(img);
dxa = [-sqrt(2)/4 0 sqrt(2)/4; -1 0 1; -sqrt(2)/4 0 sqrt(2)/4];
dya = dxa'; % deriv?e verticale
g = fspecial('gaussian', max(1, fix(5 * sigma)), sigma); % gaussien
Ixa = conv2(imd, dxa, 'same');
Iya = conv2(imd, dya, 'same');
Ixa2 = conv2(Ixa.^2, g, 'same');
Iya2 = conv2(Iya.^2, g, 'same');
Ixya = conv2(Ixa.*Iya, g, 'same');
R = Ixa2 .* Iya2 - Ixya.^2 - lambda * (Ixa2 + Iya2).^2;

R1 = (1000 / (max(max(R)))) * R; %normalisation
[u, v] = find(R1 <= seuil);
nb = length(u);

for k = 1:nb
    R1(u(k), v(k)) = 0;
end

R11 = zeros(m + 2 * r, n + 2 * r);
R11(r + 1:m + r, r + 1:n + r) = R1;

[m1, n1] = size(R11);
for i = r + 1:m1 - r
    for j = r + 1:n1 - r
        fenetre = R11(i - r:i + r, j - r:j + r);
        ma = max(max(fenetre));
        if fenetre(r + 1, r + 1) < ma
            R11(i, j) = 0;
        end
    end
end

axes(handles.axes2);
imshow(img);
hold on;
R11 = R11(r + 1:m + r, r + 1:n + r);
[x, y] = find(R11);
nb = length(x);
plot(y, x, '.r')
title('Detection des points d''interet');
% --------------------------------------------------------------------
function Erosion_Callback(hObject, eventdata, handles)
% hObject    handle to Erosion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Charger une image en niveaux de gris
    image = handles.courant_data;

    % Convertir l'image en niveaux de gris si elle est en couleur
    if size(image, 3) == 3
        image = rgb2gray(image);
    end

    % D?finir l'?l?ment structurant pour l'?rosion
    se = strel('square', 3);

    % Effectuer l'?rosion de l'image
    image_eroded = imerode(image, se);

    % Afficher l'image ?rod?e
    axes(handles.axes2);
    imshow(image_eroded);
    title('Image erodee');

% --------------------------------------------------------------------
function Dilatation_Callback(hObject, eventdata, handles)
% hObject    handle to Dilatation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Charger une image en niveaux de gris
    image = handles.courant_data;

    % Convertir l'image en niveaux de gris si elle est en couleur
    if size(image, 3) == 3
        image = rgb2gray(image);
    end

    % D?finir l'?l?ment structurant pour la dilatation
    se = strel('square', 3);

    % Effectuer la dilatation de l'image
    image_dilated = imdilate(image, se);

    % Afficher l'image dilat?e
    axes(handles.axes2);
    imshow(image_dilated);
    title('Image dilatee');

% --------------------------------------------------------------------
function Ouverture_Callback(hObject, eventdata, handles)
% hObject    handle to Ouverture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Charger une image en niveaux de gris
    image = handles.courant_data;

    % Convertir l'image en niveaux de gris si elle est en couleur
    if size(image, 3) == 3
        image = rgb2gray(image);
    end

    % Convertir l'image en type double pour les op?rations de traitement d'images
    image = double(image);

    % D?finir l'?l?ment structurant pour l'ouverture (par exemple, une grille de 5x5 avec des ?l?ments structurants carr?s)
    se = strel('square', 5);

    % Effectuer l'ouverture de l'image
    image_opened = imopen(image, se);

    % Convertir l'image de retour en uint8 pour l'affichage
    image_opened = uint8(image_opened);

    % Afficher l'image ouverte
    axes(handles.axes2);
    imshow(image_opened);
    title('Image ouverte');
% --------------------------------------------------------------------
function Fermeture_Callback(hObject, eventdata, handles)
% hObject    handle to Fermeture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 % Charger une image en niveaux de gris
    image = handles.courant_data;

    % Convertir l'image en niveaux de gris si elle est en couleur
    if size(image, 3) == 3
        image = rgb2gray(image);
    end

    % Convertir l'image en type double pour les op?rations de traitement d'images
    image = double(image);

    % D?finir l'?l?ment structurant pour la fermeture (par exemple, une grille de 5x5 avec des ?l?ments structurants carr?s)
    se = strel('square', 5);

    % Effectuer la fermeture de l'image
    image_closed = imclose(image, se);

    % Convertir l'image de retour en uint8 pour l'affichage
    image_closed = uint8(image_closed);

    % Afficher l'image ferm?e
    axes(handles.axes2);
    imshow(image_closed);
    title('Image fermee');
% -------------------------------------------------------------------
function FPB_Callback(hObject, eventdata, handles)
% hObject    handle to FPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.courant_data;
F=fftshift(fft2(I));
% %calcul de la taille de l'image;
M=size(F,1);
N=size(F,2);
P=size(F,3);
H0=zeros(M,N);
D0=1;
M2=round(M/2);
N2=round(N/2);
H0(M2-D0:M2+D0,N2-D0:N2+D0)=1;
for i=1:M
    for j=1:N
        G(i,j)=F(i,j)*H0(i,j);
    end
end
g=ifft2(G);
imshow(abs(g),[0,255])



% --------------------------------------------------------------------
function FPBB_Callback(hObject, eventdata, handles)
% hObject    handle to FPBB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% R?cup?rer l'image
I = handles.courant_data;
%I = imread('eight.tif');
F=fftshift(fft2(I));
%calcul de la taille de l'image;
M=size(F,1);
N=size(F,2);
P=size(F,3);
H0=zeros(M,N);
D0=1;
M2=round(M/2);
N2=round(N/2);
H0(M2-D0:M2+D0,N2-D0:N2+D0)=1;
n=3;
for i=1:M
    for j=1:N
        H(i,j)=1/(1+(H0(i,j)/D0)^(2*n));
        G(i,j)=F(i,j)*H0(i,j);
    end
end
g=ifft2(G);
imshow(abs(g),[0,255]);%title('image filtr?e');

% --------------------------------------------------------------------
function FPH_Callback(hObject, eventdata, handles)
% hObject    handle to FPH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% R?cup?rer l'image
I = handles.courant_data;

% Effectuer le traitement
F = fftshift(fft2(I));
[M, N] = size(F);
D0 = 1;
M2 = round(M/2);
N2 = round(N/2);

% Cr?er le filtre passe-haut
H1 = ones(M,N);
H1(M2-D0:M2+D0,N2-D0:N2+D0) = 0;

% Appliquer le filtre dans le domaine de fr?quence
G = F .* H1;

% Obtenir l'image filtr?e en revenant dans le domaine spatial
g = ifft2(ifftshift(G));

% Afficher l'image filtr?e
axes(handles.axes2);
imshow(255 - abs(g), [0, 255]);
% --------------------------------------------------------------------
function FPHB_Callback(hObject, eventdata, handles)
% hObject    handle to FPHB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% R?cup?rer l'image
I = handles.courant_data;

% Effectuer la transformation de Fourier
F = fftshift(fft2(I));
[M, N] = size(F);
D0 = 1;
M2 = round(M/2);
N2 = round(N/2);

% Cr?er le filtre passe-haut
H1 = ones(M,N);
H1(M2-D0:M2+D0,N2-D0:N2+D0) = 0;

% Appliquer le filtre dans le domaine de fr?quence
G = F .* H1;

% Appliquer la fonction de transfert
n = 3;
H = 1 ./ (1 + (H1 / D0).^(2*n));

% Appliquer le filtre modifi?
G = G .* H;

% Obtenir l'image filtr?e en revenant dans le domaine spatial
g = ifft2(ifftshift(G));

% Afficher l'image filtr?e
axes(handles.axes2);
imshow(255 - abs(g), [0, 255]);


% --------------------------------------------------------------------
function hough_Callback(hObject, eventdata, handles)
% hObject    handle to hough (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function droite_Callback(hObject, eventdata, handles)
% hObject    handle to droite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageIn=handles.courant_data;
if length( size(imageIn,3))>=3
    imageIn=rgb2gray(imageIn);
end
BW = edge(imageIn,'canny');
[H,theta,rho] = hough(BW);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW,theta,rho,P,'FillGap',5,'MinLength',7);
axes(handles.axes2);
imshow(imageIn), hold on
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    % Plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red')

% --------------------------------------------------------------------
function cercle_Callback(hObject, eventdata, handles)
% hObject    handle to cercle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 imageIn=handles.courant_data;
if length( size(imageIn,3))>=3
    imageIn=rgb2gray(imageIn);
end
e = edge(imageIn, 'canny');
radii = 15:1:40;
h = circle_hough(e, radii, 'same', 'normalise');
peaks = circle_houghpeaks(h, radii, 'nhoodxy', 15, 'nhoodr', 21, 'npeaks', 10);
axes(handles.axes2);
imshow(imageIn), hold on
for peak = peaks
    [x, y] = circlepoints(peak(3));
    plot(x+peak(1), y+peak(2), 'g-');
end
hold off


    
    
% --------------------------------------------------------------------
function CI_Callback(hObject, eventdata, handles)
% hObject    handle to CI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=handles.courant_data;
%se = strel('line',11,90);
se = strel('disk',4);
erodedI = imerode(image,se);
image=double(image)-double(erodedI);
nv=uint8(image);
axes(handles.axes2);
subimage(nv);
% --------------------------------------------------------------------
function CE_Callback(hObject, eventdata, handles)
% hObject    handle to CE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%se = strel('line',11,90);
se = strel('disk',4);
dilatedI = imdilate(image,se);
image=double(dilatedI)-double(image);
nv=uint8(image);
axes(handles.axes2);
subimage(nv);


% --------------------------------------------------------------------
function Gm_Callback(hObject, eventdata, handles)
% hObject    handle to Gm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=handles.courant_data;
%se = strel('line',11,90);
se = strel('disk',4);
erodedI = imerode(image,se);
dilatedI = imdilate(image,se);
image=double(dilatedI)-double(erodedI);
nv=uint8(image);
axes(handles.axes2);
subimage(nv);
% --------------------------------------------------------------------
function BTH_Callback(hObject, eventdata, handles)
% hObject    handle to BTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=handles.courant_data;
%element structurant de type disc avec rayon = 4 pixel
se = strel('disk',4);
F = imclose(image,se);
image=double(F)-double(image);
nv=uint8(image);
axes(handles.axes2);
subimage(nv);
% --------------------------------------------------------------------
function WTH_Callback(hObject, eventdata, handles)
% hObject    handle to WTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%figure, imhist(image);
%element structurant de type disc avec rayon = 4 pixel
se = strel('disk',4);
O = imopen(image,se);
image=double(image)-double(O);
nv=uint8(image);
axes(handles.axes2);
subimage(nv);


% --------------------------------------------------------------------
function passBande_Callback(hObject, eventdata, handles)
% hObject    handle to passBande (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
if size(image, 3) == 3
     image = rgb2gray(image);
end

F = fftshift(fft2(image));
[M, N] = size(F);
% D?finir les param?tres du filtre passe-bande
D0 = 50; % Fr?quence de coupure basse
D1 = 150; % Fr?quence de coupure haute
% Cr?er un masque pour le filtre passe-bande
H = zeros(M, N);
for u = 1:M
     for v = 1:N
        D = sqrt((u - M/2)^2 + (v - N/2)^2);
        if D >= D0 && D <= D1
             H(u, v) = 1;
        end
     end
end
G = F .* H;
g = ifft2(ifftshift(G));
axes(handles.axes2);
imshow(abs(g), []);

% --------------------------------------------------------------------
function sobel_Callback(hObject, eventdata, handles)
% hObject    handle to sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Read the image data from handles
    image = handles.courant_data;

    % Convert image to double
    image = double(image);

    % Size of the image
    [n, m] = size(image);

    % Initialize output matrices
    output = zeros(size(image));
    outputhor = zeros(size(image));
    outputver = zeros(size(image));

    % Sobel masks
    maskhor = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
    maskver = [-1, -2, -1; 0, 0, 0; 1, 2, 1];

    % Apply Sobel operator
    for i = 4:(m-3)
        for j = 4:(n-3)
            for k = 1:3
                for l = 1:3
                    outputhor(i,j) = outputhor(i,j) + image(i-k,j-l)*maskhor(k,l);
                    outputver(i,j) = outputver(i,j) + image(i-k,j-l)*maskver(k,l);
                end
            end
        end
    end

    % Compute gradient magnitude
    output = sqrt(outputhor.^2 + outputver.^2);

    % Convert the filtered image back to uint8
    output = uint8(output);

    % Display the filtered image on axes2
    axes(handles.axes2);
    subimage(output);

% --------------------------------------------------------------------
function Kirsch_Callback(hObject, eventdata, handles)
% hObject    handle to Kirsch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imageIn=handles.courant_data;
x=double(imageIn);
g1=[5,5,5; -3,0,-3; -3,-3,-3];
g2=[5,5,-3; 5,0,-3; -3,-3,-3];
g3=[5,-3,-3; 5,0,-3; 5,-3,-3];
g4=[-3,-3,-3; 5,0,-3; 5,5,-3];
g5=[-3,-3,-3; -3,0,-3; 5,5,5];
g6=[-3,-3,-3; -3,0,5;-3,5,5];
g7=[-3,-3,5; -3,0,5;-3,-3,5];
g8=[-3,5,5; -3,0,5;-3,-3,-3];
x1=imfilter(x,g1,'replicate');
x2=imfilter(x,g2,'replicate');
x3=imfilter(x,g3,'replicate');
x4=imfilter(x,g4,'replicate');
x5=imfilter(x,g5,'replicate');
x6=imfilter(x,g6,'replicate');
x7=imfilter(x,g7,'replicate');
x8=imfilter(x,g8,'replicate');
y1=max(x1,x2);
y2=max(y1,x3);
y3=max(y2,x4);
y4=max(y3,x5);
y5=max(y4,x6);
y6=max(y5,x7);
y7=max(y6,x8);
y=y7;
axes(handles.axes2);
subimage(uint8(y));


% --------------------------------------------------------------------
function marH_Callback(hObject, eventdata, handles)
% hObject    handle to marH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageIn=handles.courant_data;
im=im2double(imageIn);
%smoothening the image with a filter
gfilter= [0 0 1 0 0;
0 1 2 1 0;
1 2 -16 2 1;
0 1 2 1 0;
0 0 1 0 0];
smim=conv2(im,gfilter);
% finding the zero crossings
[rr,cc]=size(smim);
zc=zeros([rr,cc]);
for i=2:rr-1
    for j=2:cc-1
        if (smim(i,j)>0)
            if (smim(i,j+1)>=0 && smim(i,j-1)<0) || (smim(i,j+1)<0 && smim(i,j-1)>=0)
                zc(i,j)= smim(i,j+1);
            elseif (smim(i+1,j)>=0 && smim(i-1,j)<0) || (smim(i+1,j)<0 && smim(i-1,j)>=0)
                zc(i,j)= smim(i,j+1);
            elseif (smim(i+1,j+1)>=0 && smim(i-1,j-1)<0) || (smim(i+1,j+1)<0 && smim(i-1,j-1)>=0)
                zc(i,j)= smim(i,j+1);
            elseif (smim(i-1,j+1)>=0 && smim(i+1,j-1)<0) || (smim(i-1,j+1)<0 && smim(i+1,j-1)>=0)
                zc(i,j)=smim(i,j+1);
            end
        end
    end
end
otpt=im2uint8(zc);
% thresholding
otptth= otpt>105;
figure;
subplot(2,2,1);imshow(im);title('Origional image');
subplot(2,2,2);imshow(smim);title('Smoothened image');
subplot(2,2,3);imshow(otpt);title('Output image');
subplot(2,2,4);imshow(otptth);title('Output image with threshold');
% final result
figure, imshow(otptth);
axes(handles.axes2);
subimage(smim)
% --------------------------------------------------------------------
function canny_Callback(hObject, eventdata, handles)
% hObject    handle to canny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img=handles.courant_data;

T_Low = 0.075;
T_High = 0.175;
%Gaussian Filter Coefficient
B = [2, 4, 5, 4, 2; 4, 9, 12, 9, 4;5, 12, 15, 12, 5;4, 9, 12, 9, 4;2, 4, 5, 4, 2 ];
B = 1/159.* B;
%Convolution of image by Gaussian Coefficient
A=conv2(img, B, 'same');
%Filter for horizontal and vertical direction
KGx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
KGy = [1, 2, 1; 0, 0, 0; -1, -2, -1];
%Convolution by image by horizontal and vertical filter
Filtered_X = conv2(A, KGx, 'same');
Filtered_Y = conv2(A, KGy, 'same');
%Calculate directions/orientations
arah = atan2 (Filtered_Y, Filtered_X);
arah = arah*180/pi;
pan=size(A,1);
leb=size(A,2);
%Adjustment for negative directions, making all directions positive
for i=1:pan
    for j=1:leb
        if (arah(i,j)<0)
            arah(i,j)=360+arah(i,j);
        end;
    end;
end;
arah2=zeros(pan, leb);
%Adjusting directions to nearest 0, 45, 90, or 135 degree
for i = 1 : pan
    for j = 1 : leb
        if ((arah(i, j) >= 0 ) && (arah(i, j) < 22.5) || (arah(i, j) >= 157.5) && (arah(i, j) < 202.5) || (arah(i, j) >=337.5) && (arah(i, j) <= 360))
            arah2(i, j) = 0;
        elseif ((arah(i, j) >= 22.5) && (arah(i, j) < 67.5) || (arah(i, j) >= 202.5) && (arah(i, j) < 247.5))
            arah2(i, j) = 45;
        elseif ((arah(i, j) >= 67.5 && arah(i, j) < 112.5) || (arah(i, j) >= 247.5 && arah(i, j) < 292.5))
             arah2(i, j) = 90;
        elseif ((arah(i, j) >= 112.5 && arah(i, j) < 157.5) || (arah(i, j) >= 292.5 && arah(i, j) < 337.5))
             arah2(i, j) = 135;
        end;
    end;
end;
figure, imagesc(arah2); colorbar;
%Calculate magnitude
magnitude = (Filtered_X.^2) + (Filtered_Y.^2);
magnitude2 = sqrt(magnitude);
BW = zeros (pan, leb);
%Non-Maximum Supression
for i=2:pan-1
    for j=2:leb-1
        if (arah2(i,j)==0)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i,j+1), magnitude2(i,j-1)]));
        elseif (arah2(i,j)==45)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j-1), magnitude2(i-1,j+1)]));
        elseif (arah2(i,j)==90)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j), magnitude2(i-1,j)]));
        elseif (arah2(i,j)==135)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j+1), magnitude2(i-1,j-1)]));
        end;
    end;
end;
BW = BW.*magnitude2;
figure, imshow(BW);
%Hysteresis Thresholding
T_Low = T_Low * max(max(BW));
T_High = T_High * max(max(BW));
T_res = zeros (pan, leb);
for i = 1 : pan
    for j = 1 : leb
        if (BW(i, j) < T_Low)
            T_res(i, j) = 0;
        elseif (BW(i, j) > T_High)
            T_res(i, j) = 1;
%Using 8-connected components
        elseif ( BW(i+1,j)>T_High || BW(i-1,j)>T_High || BW(i,j+1)>T_High || BW(i,j-1)>T_High || BW(i-1,j-1)>T_High || BW(i-1, j+1)>T_High || BW(i+1, j+1)>T_High || BW(i+1, j-1)>T_High)
            T_res(i,j) = 1;
        end;

    end;
end;
edge_final = uint8(T_res.*255);
axes(handles.axes2);
subimage(edge_final);

% --------------------------------------------------------------------
function Gaussien55_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussien55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 % Read the image data from handles
    image = handles.courant_data;

    % Convert image to double
    image = double(image);

    % Size of the image
    [n, m] = size(image);

    % Define the Gaussian filter
    H = (1/256) * [1 4 6 4 1 ; 4 16 24 16 4 ; 6 24 36 24 6 ; 4 16 24 16 4 ; 1 4 6 4 1];

    % Apply the Gaussian filter to the image
    b = zeros(n, m);
    for x = 3 : n-2
        for y = 3 : m-2
            f = image(x-2:x+2, y-2:y+2);
            v = f .* H;
            b(x, y) = sum(v(:));
        end
    end

    % Convert the filtered image back to uint8
    b = uint8(b);

    % Display the filtered image on axes2
    axes(handles.axes2);
    imshow(b);
    title('Gaussian (5,5)');
% --------------------------------------------------------------------
function moy55_Callback(hObject, eventdata, handles)
% hObject    handle to moy55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Read the image data from handles
    image = handles.courant_data;

    % Convert image to double
    image = double(image);

    % Size of the image
    [n, m] = size(image);

    % Define the mean filter
    H = (1/25) * ones(5, 5);

    % Apply the mean filter to the image
    b = zeros(n, m);
    for x = 3 : n-2
        for y = 3 : m-2
            f = image(x-2:x+2, y-2:y+2);
            v = f .* H;
            b(x, y) = sum(v(:));
        end
    end

    % Convert the filtered image back to uint8
    b = uint8(b);

    % Display the filtered image on axes2
    axes(handles.axes2);
    imshow(b);
    title('Moyanneur (5,5)');

% --------------------------------------------------------------------
function pyramidale_Callback(hObject, eventdata, handles)
% hObject    handle to pyramidale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Read the image data from handles
    image = handles.courant_data;

    % Convert image to double
    image = double(image);

    % Size of the image
    [n, m] = size(image);

    % Define the Pyramidal filter
    H = (1/81) * [1 2 3 2 1 ; 2 4 6 4 2 ; 3 6 9 6 3 ; 2 4 6 4 2 ; 1 2 3 2 1];

    % Apply the Pyramidal filter to the image
    b = zeros(n, m);
    for x = 3 : n-2
        for y = 3 : m-2
            f = image(x-2:x+2, y-2:y+2);
            v = f .* H;
            b(x, y) = sum(v(:));
        end
    end

    % Convert the filtered image back to uint8
    b = uint8(b);

    % Display the filtered image on axes2
    axes(handles.axes2);
    imshow(b);
    title('Pyramidal Filtered Image');

% --------------------------------------------------------------------
function Decalage_Callback(hObject, eventdata, handles)
% hObject    handle to Decalage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
    [n, m] = size(image);
    image = double(image);
    output = image;
    L = 100;
    
    % Apply histogram stretching
    [l, c] = size(image);
    v = image;
    for i = 1:l
        for j = 1:c
            fpixel = image(i, j) + L;
            % Check if pixel value is within [0, 255]
            if fpixel > 255
                fpixel = 255;
            elseif fpixel < 0
                fpixel = 0;
            end
            v(i, j) = fpixel;
        end
    end
    
    % Display stretched image on axes2
    axes(handles.axes2);
    imshow(uint8(v), [0, 255]);

% --------------------------------------------------------------------
function multiple_Callback(hObject, eventdata, handles)
% hObject    handle to multiple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image = handles.courant_data;
    
    % Adjust brightness and contrast
    P = 1.5;
    adjusted_image = double(image) * P;
    adjusted_image(adjusted_image > 255) = 255;
    adjusted_image(adjusted_image < 0) = 0;
    
    % Display adjusted image on axes2
    axes(handles.axes2);
    imshow(uint8(adjusted_image), [0, 255]);
    title('Adjusted Image');
    axis off;
% --------------------------------------------------------------------
function seuillage_Callback(hObject, eventdata, handles)
% hObject    handle to seuillage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageIn=handles.courant_data;
thresholded = thresholdLocally(imageIn);
axes(handles.axes2);
imshow(thresholded);

% --------------------------------------------------------------------
function Prewitt_Callback(hObject, eventdata, handles)
% hObject    handle to Prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
output=image;
%image=rgb2gray(image);
output=zeros(size(image));
outputhor=zeros(size(image));
outputver=zeros(size(image));
maskhor = [-1,0,1;-1,0,1;-1,0,1];
maskver = [-1,-1,-1;0,0,0;1,1,1];
for i=4:(m-3)
    for j=4:(n-3)
        for k=1:3
            for l=1:3
                outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);
                outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);
            end
        end
    end
end
for i=4:(m-3)
    for j=4:(n-3)
        output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j));
    end
end
output=uint8(output);
axes(handles.axes2);
subimage(output)