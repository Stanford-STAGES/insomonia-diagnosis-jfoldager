function out = hjort_parameters(sig,fs,segmentLength)
%%  Calculates Hjort's 3 parameters columnwise from sig divided into segmentLength 
    x = sig(:);
    x         = reshape(x(1:end-mod(length(x) ...
                ,fs*segmentLength)),fs*segmentLength,[]);
    out.act = activity(x);
    dt      = repmat(1/fs,size(x,1)-1,1);
    for i = 1:size(x,2)
        dxV         = diff([x(:,i)]);
        dxdt        = dxV./dt;
%         plot(x(:,i)/max(abs(x(:,i)))),hold on,plot(dxV/max(abs(dxV))),hold on,plot(dxdt/max(abs(dxdt)));
%         legend('sig','dsig','dsig/dt')
        out.mob(i)  = sqrt(activity(dxdt)/out.act(i));
        ddxV        = diff([0;dxV]);
        mx2         = mean(x(:,i).^2);
        mdx2        = mean(dxV.^2);
        mddx2       = mean(ddxV.^2);
        mob         = mdx2 / mx2;
        out.mob(i)  = sqrt(mddx2 / mdx2 - mob);
        out.com(i)  = sqrt(mob);
    end
    out.ALL = [out.act',out.mob',out.com'];
end
function act = activity(x)
    act = var(x,0,1);
end