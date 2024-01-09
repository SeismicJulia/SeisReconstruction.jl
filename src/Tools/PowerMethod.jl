function PowerMethod(x0,operators,parameters)
    
    x= x0;
    α=0.0;
    for k = 1:10;
        aux=LinearOperator(x,operators,parameters,adj=false)
        y=LinearOperator(aux,operators,parameters,adj=true)
        n = norm(y,2);
        x = y/n;
        α = n;
    end
    return α
end
