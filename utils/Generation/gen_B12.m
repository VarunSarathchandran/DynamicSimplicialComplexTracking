function out = gen_B12(N)
    NN = 1:N;
    A = ones(N)-eye(N);
    E = N*(N-1)/2;
    T = trace(A^3)/6;
    G = graph(A);
    B1 = full(incidence(G));%incidence matrix
    B2 = zeros(E,T);
    t = 1;
    e = 1;
    for i = NN
        for j = NN(NN>i)
            for k = NN(NN>j)
                idx1 = N*(i-1)+j-sum(1:i);
                idx2 = N*(i-1)+k-sum(1:i);
                idx3 = N*(j-1)+k-sum(1:j);
                %disp(['Triangulo: ' num2str(i) num2str(j) num2str(k)])
                B2([idx1,idx3],t) = 1;
                B2([idx2],t) = -1;
                t = t+1;
            end
            e = e+1;
        end
    end
    out.B1 = B1;
    out.B2 = B2;
end