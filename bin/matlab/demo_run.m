function demo_run
fbase_filename = '~/data/meshtk_workshop/sample.fbase';
eigenvalue_filename = '~/data/meshtk_workshop/sample.ev/_ev.ascii';

[L D] = TriMeshTKFBaseRead(fbase_filename, eigenvalue_filename);
L = L(1:100,1:100);
D = D(1:101);

tmp = log(3*D(end)/D(2));
logtime_step = tmp/99;
T=1/D(end)*exp(0:logtime_step:tmp);

time_c = length(T);


numberofeigenvalues = 100;
E=zeros(numberofeigenvalues, time_c);

evtrack = zeros(numberofeigenvalues, time_c);

numberoftracks = 100;

for t = time_c:-1:1
    
    %f = @(x) exp(-x*T(t));
    f=@(x) x*T(t).*exp(-x*T(t));
    df = @(x) (1-x*T(t)) .* x .* exp(-x*T(t));
    
    
    
    
    
    if (t==time_c)
        M = TriMeshTKFBase2M(f, L, D);
        dM = TriMeshTKFBase2M(df, L, D);
        [curEigVector, curEigValue] = eigs(M, numberofeigenvalues);
        U0 = (curEigVector(:,1)' * dM * curEigVector(:,1))/abs(curEigValue(1,1));
        
        
        E(:,t) = diag(curEigValue);        
        E(:,t)=log(abs(E(:,t)/E(1,t)));
        %evtrack(1,t) =1;
        evtrack(end-numberoftracks+1:end,t) = 1:numberoftracks;
    else

        preEigVector = curEigVector;
        preEigValue = curEigValue;
        V0 =U0;
        predM = dM;
        
        
        M = TriMeshTKFBase2M(f, L, D);
        dM = TriMeshTKFBase2M(df, L, D);
        [curEigVector, curEigValue] = eigs(M, numberofeigenvalues);
        U0 = (curEigVector(:,1)' * dM * curEigVector(:,1))/abs(curEigValue(1,1));
        E(:,t) = diag(curEigValue);
        
        E(:,t)=log(abs(E(:,t)/E(1,t)));
%         for j=1:numberofeigenvalues
% 
%             if (evtrack(j, t+1) ~= 0)
%                 index = track_eigvector(j);
%                 if (index ~=0)
%                     if (evtrack(j,t+1) > evtrack(index,t))
%                         evtrack(index, t) = evtrack(j, t+1);
%                     end
%                 end
%             end
%         end
        score = - ones(numberofeigenvalues,numberoftracks);
        for j=1:numberoftracks
            tmp = 1:numberofeigenvalues;
            node = tmp(evtrack(:,t+1) == j);
            if (length(node)==1)
                track_eigvector(node,j);
            end
        end
        
        for j=1:numberoftracks
            [tmp, tmpindex] = max(score);
            [energy, index] = max(tmp);
            if (energy > 0)
                evtrack(tmpindex(index),t) = index;
            else
                break;
            end
            score(:, index) = -1;
            score(tmpindex(index),:) = -1;
        end
    end
    
    
    
end

logT=log(T);
%logT = T;

count = ones(numberofeigenvalues, 1) * (1:numberofeigenvalues); 

%for i=numberoftracks-20:numberoftracks
%    curv = E(evtrack == i);
%    pT = logT(count(evtrack ==i));
%    plot(pT,curv,'-ob','MarkerSize',2);
%    
%    hold on;
%end

for i=1:numberofeigenvalues
    plot(logT,E(i,:),'o', 'MarkerSize',2);
    hold on;
end
ylim([-log(numberofeigenvalues) 0]);
xlim([min(logT),max(logT)]);

    
    function  track_eigvector(k,classno)        
        V = preEigVector(:,k);
        Vev = abs(preEigValue(k,k));        
        dV = T(t+1)* ( V'*predM*V/Vev - V0);
        Vev = E(k,t+1);%log(Vev / preEigValue(1,1));
                
        for ii=1:numberofeigenvalues
            U = curEigVector(:,ii);
            Uev = abs(curEigValue(ii,ii));
            dU = T(t)*( U'*dM*U/Uev - U0);
            Uev = E(ii,t);%log(Uev / curEigValue(1,1));
            
            vecsim = 2*abs(acos(abs(U'*V)))/pi;
            gradsim = 3*abs(atan((dU-dV)/(1+dU*dV)))/pi;
            disgrad = (Vev - Uev)/logtime_step;
            evsim = abs(atan((disgrad -(dV+dU)/2)/(1+disgrad*(dV+dU)/2)))/pi;            
            
            score(ii,classno) = 1- vecsim - gradsim - evsim;
        end
    end
end