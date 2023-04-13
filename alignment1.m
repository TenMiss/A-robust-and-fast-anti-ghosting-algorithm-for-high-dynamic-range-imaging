function [delta, phi] = alignment1(im)
im0{1} = double(im{1});
lp = fspecial('ga',5,1);
for i=2:3
    im0{i} = imresize(conv2(double(im0{i-1}),lp,'same'),0.5,'bilinear');
end
im1{1} = double(im{2});
for i=2:3
    im1{i} = imresize(conv2(double(im1{i-1}),lp,'same'),0.5,'bilinear');
end
stot = zeros(1,3);
for pyrlevel=3:-1:1
    f0 = im0{pyrlevel};
    f1 = im1{pyrlevel};
    [y0,x0]=size(f0);
    xmean=x0/2; ymean=y0/2;
    x=kron([-xmean:xmean-1],ones(y0,1));
    y=kron(ones(1,x0),[-ymean:ymean-1]');
    sigma=1;
    g3 = exp(-(y.^2+x.^2)/(2*sigma^2))/2/pi/sigma^2;
    g1 = -g3.*y; 
    g2 = -g3.*x; 
    a=real(ifft2(fft2(f1).*fft2(g2))); 
    c=real(ifft2(fft2(f1).*fft2(g1))); 
    b=real(ifft2(fft2(f1).*fft2(g3)))-real(ifft2(fft2(f0).*fft2(g3))); 
    R=c.*x-a.*y; 
    a11 = sum(sum(a.*a)); a12 = sum(sum(a.*c)); a13 = sum(sum(R.*a));
    a21 = sum(sum(a.*c)); a22 = sum(sum(c.*c)); a23 = sum(sum(R.*c)); 
    a31 = sum(sum(R.*a)); a32 = sum(sum(R.*c)); a33 = sum(sum(R.*R));
    b1 = sum(sum(a.*b)); b2 = sum(sum(c.*b)); b3 = sum(sum(R.*b));
    Ainv = [a11 a12 a13; a21 a22 a23; a31 a32 a33]^(-1);
    s = Ainv*[b1; b2; b3];
    st = s;
    it=1;
    while ((abs(s(1))+abs(s(2))+abs(s(3))*180/pi/20>0.1)&it<25)
        f0_ = shift(f0,-st(1),-st(2));
        f0_ = imrotate(f0_,-st(3)*180/pi,'bilinear','crop');
        b = real(ifft2(fft2(f1).*fft2(g3)))-real(ifft2(fft2(f0_).*fft2(g3)));
        s = Ainv*[sum(sum(a.*b)); sum(sum(c.*b)); sum(sum(R.*b))];
        st = st+s;
        it = it+1;
    end
    st(3)=-st(3)*180/pi;
    st = st';
    st(1:2) = st(2:-1:1);
    stot = [2*stot(1:2)+st(1:2) stot(3)+st(3)];
    if pyrlevel>1
        im1{pyrlevel-1} = imrotate(im1{pyrlevel-1},-stot(3),'bilinear','crop');
        im1{pyrlevel-1} = shift(im1{pyrlevel-1},2*stot(2),2*stot(1));             
    end
end
phi = stot(3);
delta = stot(1:2);
