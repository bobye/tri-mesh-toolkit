fbase_filename = '~/data/meshtk_workshop/T77.fbase';
eigenvalue_filename = '~/data/meshtk_workshop/T77.ev/_ev.ascii';

[L D] = TriMeshTKFBaseRead(fbase_filename, eigenvalue_filename);

tmp = log(4*D(end)/D(2));
T=1/D(end)*exp(0:tmp/99:tmp);
time_c = length(T);
E=zeros(length(T),100);

for t = 1:time_c
    %f = @(x) exp(-x*T(t));
    f=@(x) x*T(t).*exp(-x*T(t));
    M = TriMeshTKFBase2M(f,L,D);
    E(t,:) = eigs(M,100); E(t,:)= log(abs(E(t,:)/E(t,1)));
end

logT=1:time_c;

%plot(logT,E(:,2),logT,E(:,3),logT,E(:,4),logT,E(:,5),logT,E(:,6),logT,E(:,7),logT,E(:,8),logT,E(:,9),logT,E(:,10),logT,E(:,11),logT,E(:,12),logT,E(:,13),logT,E(:,14),logT,E(:,15),logT,E(:,16),logT,E(:,17),logT,E(:,18),logT,E(:,19),logT,E(:,20));
for i=2:100
    plot(logT,E(:,i),'o','MarkerSize',2);
    hold on;
end
ylim([-log(length(D)) 0]);