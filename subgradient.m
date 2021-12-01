function [realPP] = subgradient(f)

%% set up
N = length(f);
%% connectivity
epsi = 1e-16;
CC.NumObjects = 1;
CC.PixelIdxList{1}=1;
alternativeSignal = f(1);
for iii=2:1:N
    if(abs(f(iii)-f(iii-1))<epsi)
        CC.PixelIdxList{end}=[CC.PixelIdxList{end},iii];
    else
        alternativeSignal(end+1)=f(iii);
        CC.NumObjects = CC.NumObjects + 1;
        CC.PixelIdxList{end+1} = iii;
    end
end
gradiLyRight = [sign(alternativeSignal(2:1:end)-alternativeSignal(1:1:end-1)),0];
gradiLyLeft = [0,sign(alternativeSignal(2:1:end)-alternativeSignal(1:1:end-1))];
pp = -gradiLyRight+gradiLyLeft;

ind = 1;
for iii=1:1:CC.NumObjects
    for jjj=1:1:length(CC.PixelIdxList{iii})
        realPP(ind)=pp(iii)/length(CC.PixelIdxList{iii});
        ind = ind + 1;
    end
end
return;
end
