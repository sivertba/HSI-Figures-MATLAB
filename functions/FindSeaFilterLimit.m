function limit = FindSeaFilterLimit(landEst)

[N,edges]=histcounts(landEst);

gv=[];
gvCount=0;

startIndex=find(N==max(N(edges<0)));

for i = startIndex:(length(N)-1)
    gv(i)=N(i+1)-N(i);

    if (gv(i) >= 0) || (abs(gv(i))/N(i) <= 0.1)
        gvCount=gvCount+1;
        if gvCount>1
            break;
        end
        
    elseif gv(i)<0
        gvCount=0;
    end
    
end

limit=edges(i+1);