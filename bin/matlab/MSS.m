function ImageCode = MSS(fbase_filename, eigenvalue_filename)
%Multiscale Spectra Signature

%fbase_filename = '~/data/meshtk_workshop/mesh2.fbase';
%eigenvalue_filename = '~/data/meshtk_workshop/mesh2.ev/_ev.ascii';
[L D] = TriMeshTKFBaseRead(fbase_filename, eigenvalue_filename);

L = L(1:100,1:100);
D = D(1:101);

numberofsteps = 100;
numberofeigenvalues = 100;

tmp = log(1/D(2));
logtime_step = -log(numberofeigenvalues)/(numberofsteps-1);
T=exp(tmp:logtime_step:tmp-log(numberofeigenvalues));
time_c = length(T);



E=zeros(numberofeigenvalues, time_c);

for t = time_c:-1:1
    f=@(x) x*T(t).*exp(-x*T(t));
    M = TriMeshTKFBase2M(f, L, D);
    E(:,t) = eigs(M, numberofeigenvalues);             
    E(:,t)=log(abs(E(:,t)/E(1,t)));
end
%logT=log(T);
% 
% for i=2:numberofeigenvalues
%     plot(numberofsteps:-1:1,E(i,:),'or', 'MarkerSize',2);
%     hold on;
% end
% 
% ylim([-log(numberofeigenvalues) 0]);
threshold = log(numberofeigenvalues);
pixalheight = threshold/numberofsteps;
ImageCode = zeros(numberofsteps,numberofsteps);

for t = 1:numberofsteps
    sortE = sort((E(:,t)+threshold)/pixalheight, 'descend');    
    tag = sum(sortE>0);
    if (tag == length(sortE))
        sortE(end+1) = -sortE(end);
    end
    
    for j = 1:numberofsteps
        while (j-0.5>sortE(tag)) 
            tag = tag -1;
        end
        
        if (j-.5-sortE(tag+1) > sortE(tag) -j +.5)
            ImageCode(j,t) = sortE(tag) - j +.5;
        else
            ImageCode(j,t) = j-.5 - sortE(tag+1);
        end
        
    end
end

ImageCode = ImageCode / numberofsteps;
%imshow(1-ImageCode);

% for t = 1:time_c
%     for j=2:numberofeigenvalues
%         s = E(j,t) + threshold;
%         if (s>0) 
%             i=floor(s/pixalheight+.5);
%             if (i<numberofsteps) 
%                 ImageCode(i+1,t) = ImageCode(i+1,t) + s - i*pixalheight;
%             end
%             if (i>0)
%                 ImageCode(i,t) = ImageCode(i,t) + (i+1)*pixalheight -s;
%             end
%         end
%     end
% end
% 
% ImageCode(ImageCode>1) = 1;
% ImageCode = 1 - ImageCode;

end