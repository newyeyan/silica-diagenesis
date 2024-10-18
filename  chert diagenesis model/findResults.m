function [results] = findResults(massStruct,ratioStruct,z_real,temp_vec)

% PJF 2022: find the depth, d18O and temperature at which the conversion to quartz is complete

qtzSite = find(massStruct.quartz_Si > (0.995*max(massStruct.opalCT_Si)),1,'first');

if isempty(qtzSite)
    ind = find(massStruct.quartz_Si == max(massStruct.quartz_Si));
    depth = z_real(ind);
    d18O = ratioStruct.quartz_d18O(ind);
    d30Si = ratioStruct.quartz_d30Si(ind);
    temp = temp_vec(ind);  
    if size(ind,2) > 1
        results(1) = NaN;
        results(2) = NaN;
        results(3) = NaN;
         results(4) = NaN;
         results(5) = NaN;        
    else
        
        results(1) = depth;
        results(2) = d18O;
        results(3) = temp;
         results(4) = d30Si;
         results(5) = ind;  
    end
    
else
    ind = qtzSite;
    depth = z_real(ind);
    d18O = ratioStruct.quartz_d18O(ind);
    d30Si = ratioStruct.quartz_d30Si(ind);
    temp = temp_vec(ind);
    results(1) = depth;
    results(2) = d18O;
    results(3) = temp;
     results(4) = d30Si;
    results(5) = ind;
   
end


