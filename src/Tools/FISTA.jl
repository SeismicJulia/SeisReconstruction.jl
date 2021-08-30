"""
    FISTA(x0,y,Hop,PARAM,mu,Nit)
FISTA: Solves the l2-l1 problem via Fast Iterative Shrinkage-Thresholdng Algorithm
       Given a linear operator H and it's adjoint H', the algorithm minimizes
       J = ||H x - y||_2^2 + mu ||x||_1, where H is the  linear operator
       encapsulated in Hop
# Arguments
- `y`:data
- `Hop`:linear operator
- `PARAM`:parameters to run Hop and it's adjoint
- `mu`:trade-off parameter
- `x0`:initial sol just to get size of x
Reference:  Beck and Teboulle,  2009, A Fast Iterative Shrinkage-Thresholding Algorithm
            for Linear Inverse Problems∗ SIAM J. Imaging Science,  Vol 2 (1), 183-202

"""
function FISTA(x0,y,Hop,PARAM,mu,Nit)
    x0 = randn(size(x0));
    alpha  = 1.05*power_method(x0,Hop, PARAM);
    x = zeros(size(x0));
    T = mu/(2*alpha);
    t = 1;
    yk = copy(x);
    for k=1:Nit;
        tmpx = copy(x);
        Yx = Hop(yk,PARAM,1);
        Hx=Hop( y-Yx,PARAM,-1)/alpha
        x= softThresh(yk.+Hx, T);
        tmpt = t;
        t = (1+sqrt(1+4*t^2))/2;
        yk = x + (tmpt-1)/t*(x-tmpx);
    end
    return x
end



function softThresh(x, t)
    tmp = abs.(x) .- t
    tmp = (tmp .+ abs.(tmp)) / 2
    y   = sign.(x) .* tmp
    return y
end
