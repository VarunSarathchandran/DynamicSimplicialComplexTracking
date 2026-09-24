function tf = isConfigured(cfg)
%ISCONFIGURED True when the three measurement-noise variances are supplied.

tf = ~isempty(cfg.noiseVariances);
end
