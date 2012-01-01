fbase_file = '~/data/meshtk_workshop/sample.fbase';
ev_file = '~/data/meshtk_workshop/sample.ev/_ev.ascii';

[L D] = TriMeshTKFBaseRead(fbase_file, ev_file);
ConvN = 30;
EVN = zeros(ConvN,0) ;count =1;
for i=ConvN:5:200
   M = TriMeshTKFBase2M(@(x) 1./x.^2, L(1:i,1:i), D(1:i+1));
   EV = sort(abs(eig(M)),'descend');
   EVN(:,end+1) = EV(2:ConvN+1)/EV(1);
end

M = TriMeshTKFBase2M(@(x) 1./x.^2, L, D);
EV = sort(abs(eig(M)),'descend'); EV = EV(2:ConvN+1)/EV(1);

count = size(EVN,2);
ConvChk = zeros(count,1);
for i=1:count
    ConvChk(i) = max((EVN(:,i)-EV).^2./(EVN(:,i).*EV));
end

plot(30:5:200, log(ConvChk));