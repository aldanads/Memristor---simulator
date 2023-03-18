function [R_ratio]=statistics_R_ratio(R_ratio)
m=size(R_ratio,2);

R_ratio(5,1:m/2)=R_ratio(3,1:2:m);
R_ratio(6,1)=sum(R_ratio(3,1:2:m)>1);
R_ratio(6,2)=m/2-R_ratio(6,1);
R_ratio(6,3)=mean(R_ratio(3,1:2:m));
R_ratio(6,4)=max(R_ratio(3,1:2:m));
R_ratio(6,5)=min(R_ratio(3,1:2:m));
R_ratio(6,6)=std(R_ratio(5,1:m/2));
R_ratio(6,7)=std(abs(R_ratio(5,1:m/2-1)-R_ratio(5,2:m/2)));
R_ratio(7,1:m/2)=R_ratio(3,2:2:m);
R_ratio(8,1)=sum(R_ratio(3,2:2:m)>1);
R_ratio(8,2)=m/2-R_ratio(8,1);
R_ratio(8,3)=mean(R_ratio(3,2:2:m));
R_ratio(8,4)=max(R_ratio(3,2:2:m));
R_ratio(8,5)=min(R_ratio(3,2:2:m));
R_ratio(8,6)=std(R_ratio(7,1:m/2));
R_ratio(8,7)=std(abs(R_ratio(7,1:m/2-1)-R_ratio(7,2:m/2)));

end