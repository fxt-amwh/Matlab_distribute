function xyz=my_Lhtoxyz_all(Lh)
xyz=Lh;
for i=1:size(Lh,1)
    xyz(i,:)=my_Lhtoxyz(Lh(i,:));
end
