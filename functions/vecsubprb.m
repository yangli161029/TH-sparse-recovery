function z = vecsubprb(g,lambda,L)
[Height,Width] = size(g);             
imL = Height*Width;
A  = speye(imL) + lambda * L;
z  = A\g(:);
end