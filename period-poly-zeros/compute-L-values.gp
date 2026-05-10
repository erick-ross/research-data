default(parisize, 8*10^9);



k_ub(N) = {
    if (     N == 1,           return(76););
    if (     N == 2,           return(60););
    if (  3 <= N && N <= 4,    return(56););
    if (     N == 5,           return(44););
    if (  6 <= N && N <= 7,    return(40););
    if (     N == 16,          return(102););
    if (  8 <= N && N <= 18,   return(36););
    if ( 19 <= N && N <= 29,   return(28););
    if ( 30 <= N && N <= 39,   return(24););
    if ( 40 <= N && N <= 103,  return(20););
    if (104 <= N && N <= 128,  return(16););
}


k_lb(N) = {
    if (       N == 2,         return(2 ););
    if (  1 <= N && N <= 32,   return(4 ););
    if ( 33 <= N && N <= 103,  return(12););
    if (104 <= N && N <= 128,  return(14););
}



{
for (N=1, 128, 
    for (k2=k_lb(N)\2, k_ub(N)\2-1, 
        k = 2*k2;
        dim = mfdim([N,k],0);
        printf("{'N': %d, 'k': %d, 'dim': %d}\n", N, k, dim);
        if (dim < 1, next);
        
        \\ ########## construct the newforms
        mf = mfinit([N,k],0);
        newforms = mfeigenbasis(mf);
        for(a=1, length(newforms),
            f = newforms[a];
            LLL = lfunmf(mf, f);
            \\ If there is only one embedding, then LLL is just the L-function itself,
            \\ not an array of L-functions.
            is_Lfunc = (type(LLL[length(LLL)])=="t_INT");
            L_funcs = if(is_Lfunc, [LLL], LLL);
            for(b=1, length(L_funcs),
                L_func = L_funcs[b];
                \\ ########## Now create the L-function values
                print( vector(k-1, n, lfun(L_func,n)) );
            );
        );
    );
);
}

quit;             

