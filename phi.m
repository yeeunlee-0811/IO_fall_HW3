function result = phi(x,H,k,n) 

fun = @(t) (t.^(n-k)).*((1-t).^(k-1));

result = H-(factorial(n)/(factorial(k-1)*factorial(n-k)))*integral(fun,0,x);

end

