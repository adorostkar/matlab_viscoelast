% Modification of the advection element matrix

function A2_elem = advect_modified(A1_elem)

s = size(A1_elem,1);
threshold = 0.1;
A2_elem   = A1_elem;

%for k = 1:s,
%    for l=1:s,
%%        if (k~=l)&(abs(A2_elem(k,l))<threshold), A2_elem(k,l)=0; end
%    end
%end

for k = 1:s,
    [sort_row(k,1:s),col(k,1:s)]=sort(abs(A2_elem(k,:))); 
    max_gap(k)=sort_row(k,s)-sort_row(k,s-1);
end

for k = 1:s,
    for l=s-1:-1:1,
 if (k~=col(k,l))&...
    (abs(A2_elem(k,col(k,l))) + max_gap(k)< sort_row(k,s)),...
       A2_elem(k,col(k,l))=0;  end
    end
end

return
