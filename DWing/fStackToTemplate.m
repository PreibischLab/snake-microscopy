

function template  =   fStackToTemplate(Istack)

Istack = double(Istack);
for i=1:size(Istack,3)
    X = reshape(Istack(:,:,i),[],1);
    dev = std(X);
    m  = mean(X);
    Istack(:,:,i) = (Istack(:,:,i) -m)/dev;
end

template = mean(Istack,3);


end