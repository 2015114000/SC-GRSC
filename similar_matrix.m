function [W]=similar_matrix(img,hp)

X=img';
mX        =   sum(X.^2,2)/2;
[h,w]=size(X);
W=zeros(h,h);

for off=1:h
    idx=[off:h];
    dis= ((mX(idx, :) + mX(off, :) - X(idx, :)*X(off, :)'))*2./w; 
    dis=exp(-dis./hp);
    dis(1)=0;
    dis(1)=max(dis);
    W(off,idx)=dis;
    W(idx,off)=dis;
end

end

