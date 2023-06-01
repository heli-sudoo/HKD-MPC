function plot_contact_shedule(ctactSeq)
len = length(ctactSeq);
X = categorical({'FR','FL','HR','HL'});
X = reordercats(X,{'FR','FL','HR','HL'});
Y = repmat(ones(1,len), 4, 1);
figure
b = barh(X, Y, 'stacked','FaceColor','flat','EdgeColor','none');
for i=1:len
    for l = 1:4
        b(i).CData(l,:) = ones(1,3)-ctactSeq(i,l)*ones(1,3);
        drawnow
    end    
end

end
