function out = compute_B2(A,B,w1)
    [N,E] = size(B);


    B02 = gen_B12(N).B2;   % (E × T_all 

    edge_mask = (w1 ~= 0);                 
    B2_sel = diag(edge_mask) * B02;

    w2 = (sum(abs(B2_sel), 1) == 3);       % triangles that exist 

    out.w2  = w2;               
   
end