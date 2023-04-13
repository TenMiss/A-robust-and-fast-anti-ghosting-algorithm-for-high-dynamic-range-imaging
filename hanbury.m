function [h,s,y] = hanbury(r,g,b)
%RGB2IHLS Convert red-green-blue colors to hue-saturation-luminance.
%   H = RGB2IHLS(M) converts an RGB color map to an HSY color map.
%   Each map is a matrix with any number of rows, exactly three columns,
%   and elements in the interval 0 to 1.  The columns of the input matrix,
%   M, represent intensity of red, blue and green, respectively.  The
%   columns of the resulting output matrix, H, represent hue, saturation
%   and color value, respectively.
%
%   HSY = RGB2IHLS(RGB) converts the RGB image RGB (3-D array) to the
%   equivalent HSY image HSY (3-D array).
%
%   CLASS SUPPORT
%   -------------
%   If the input is an RGB image, it can be of class uint8, uint16, or 
%   double; the output image is of class double.  If the input is a 
%   colormap, the input and output colormaps are both of class double.
% 
%   See also HSV2RGB, COLORMAP, RGBPLOT. 

%   Undocumented syntaxes:
%   [H,S,Y] = RGB2IHLS(R,G,B) converts the RGB image R,G,B to the
%   equivalent HSY image H,S,Y.
%
%   HSY = RGB2IHLS(R,G,B) converts the RGB image R,G,B to the 
%   equivalent HSV image stored in the 3-D array (HSY).
%
%   [H,S,Y] = RGB2IHLS(RGB) converts the RGB image RGB (3-D array) to
%   the equivalent HSY image H,S,Y.

% Copyright Wu Shiqian
% 17 March 2010
% Institute for Infocomm Research

switch nargin
  case 1,
     if isa(r, 'uint8'), 
        r = double(r) / 255; 
     elseif isa(r, 'uint16')
        r = double(r) / 65535;
     end
  case 3,
     if isa(r, 'uint8'), 
        r = double(r) / 255; 
     elseif isa(r, 'uint16')
        r = double(r) / 65535;
     end
     
     if isa(g, 'uint8'), 
        g = double(g) / 255; 
     elseif isa(g, 'uint16')
        g = double(g) / 65535;
     end
     
     if isa(b, 'uint8'), 
        b = double(b) / 255; 
     elseif isa(b, 'uint16')
        b = double(b) / 65535;
     end
     
  otherwise,
      error('Wrong number of input arguments.');      
end
  
threeD = (ndims(r)==3); % Determine if input includes a 3-D array

if threeD,
  g = r(:,:,2); b = r(:,:,3); r = r(:,:,1);
  siz = size(r);
  r = r(:); g = g(:); b = b(:);
elseif nargin==1,
  g = r(:,2); b = r(:,3); r = r(:,1);
  siz = size(r);
else
  if ~isequal(size(r),size(g),size(b)), 
    error('R,G,B must all be the same size.');
  end
  siz = size(r);
  r = r(:); g = g(:); b = b(:);
end
% Here be the algorithm
%luminance
%y=0.2125*r + 0.7154*g + 0.0721*b;
y=0.299*r + 0.587*g + 0.114*b;
%y=(r + g + b)/3.0;
%y=y*255;

C1 = zeros(size(y));
C2 = zeros(size(y));
C = zeros(size(y));

%hue
C1=r - 0.5*g - 0.5*b;
C2=-sqrt(3.0)/2.0*g + sqrt(3.0)/2.0*b;

C=sqrt(C1.^2 + C2.^2);

s = zeros(size(y));
h = zeros(size(y));

indic=find(C == 0);
h(indic)=0;
% faster vectorized expressions which 
% may speed up total execution by a factor of approx 2.5
% thomas knudsen, thk@kms.dk 2003-11-24
indic = find(C ~= 0 & C2 <= 0);
h(indic) = acos(C1(indic) ./ C(indic)); 
indic=find(C ~= 0 & C2 > 0);
h(indic)=2.*pi - acos(C1(indic) ./ C(indic)); 
h = h/(2*pi);
%saturation
%k=h./(pi/3.0);
%k=floor(k);
%hstar=h - k.*(pi/3.0);
%s=(2.0*C.*sin(-hstar + (2.0*pi/3.0)))/sqrt(3.0);

% simpler saturation expression which produces the same as above
s=max(r,max(g,b)) - min(r,min(g,b));
%convert h to degrees
%h=h./pi;
%h=h.*180.0;
%here ends the algorithm
if nargout<=1,
  if (threeD | nargin==3),
    h = reshape(h,siz);
    s = reshape(s,siz);
    y = reshape(y,siz);
    h=cat(3,h,s,y);
  else
    h=[h s y];
  end
else
  h = reshape(h,siz);
  s = reshape(s,siz);
  y = reshape(y,siz);
end
