function out = generate_graph_signals0(L,prms)
    if nargin < 2
       prms.M = 1e3;
       prms.sigma = 0;
       prms.verbose = false;
       prms.sig_type = 'smooth';
       prms.norm_noise = true;
    end

    if isfield(prms,'M1'); M = prms.M0; 
    else; M = 1e2; end

    if isfield(prms,'sigma'); sigma = prms.sigma; 
    else; sigma = 0; end

    if isfield(prms,'sig_type'); sig_type = prms.sig_type; 
    else; sig_type = 'smooth';end

    if isfield(prms,'verbose');verbose = prms.verbose; 
    else; verbose = false; end
   
    if isfield(prms,'norm_noise'); norm_noise = prms.norm_noise; 
    else; norm_noise = true;end

    N = size(L,1);
    h_mu = zeros(N,1);
    
    alpha = 10;%alpha indicates how smooth the signal is

 %   switch sig_type
   %     case 'smooth'
            [V, Lambda] = eig(L); 
            Lambda_inv = inv(alpha*Lambda+eye(N));
            %Lambda_inv = diag(1./(1/alpha + diag(Lambda).^2));
            %Lambda_inv = pinv(Lambda);
            H = mvnrnd(h_mu, Lambda_inv, M)';
            X_0 = V*H;
            C = pinv(L);
   %     otherwise
    %        error('ERR: Unkown signal type')
   % end
    
    p_x = norm(X_0, 'fro')^2/M;
    
    % Set sigma for normalizing the noise power
    if norm_noise
        sigma_n = sqrt(sigma*p_x/N);
    else
        sigma_n = sigma;
    end
    noise = randn(N, M)*sigma_n;
    
    X = X_0 + noise;
    Cs = X*X'/M; 
    
%     figure(1);
%     subplot(221);imagesc(X_0);colorbar();title('X0')
%     subplot(222);imagesc(noise);colorbar();title('noise');
%     subplot(223);imagesc(C);colorbar();title('C')
%     subplot(224);imagesc(Cs);colorbar();title('Cs');

    p_n = norm(noise, 'fro')^2/M;
    snr = p_x/p_n;

    if verbose
        eq_sigma = sigma^2*N/p_x;
        disp(['Mean Np: ' num2str(p_n) '   norm sigma: ' num2str(eq_sigma)])
        disp(['SNR(nat|dB): ', num2str(snr) ' | ' num2str(10*log10(snr))])
        disp(['Mean smoothness: ' num2str(trace(X_0'*L*X_0)/M)])
    end
    
    out.X_0 = X_0;
    out.X = X;
    out.C = C;
    out.Cs = Cs;
    out.snr = snr;
    out.sigma_n = sigma_n;
    out.px = p_x;
    out.pn = p_n;
    out.Lambda = Lambda;
    out.Lambda_inv = Lambda_inv;
end