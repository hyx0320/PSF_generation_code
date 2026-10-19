function [p,err] = polyfit1_rmoutliers(x,y,ex)
tmp_p = polyfit(x,y,1);
ye1 = polyval(tmp_p,x);
err1 = y - ye1;
x_tmp1 = x;
y_tmp1 = y;
x_tmp1(abs(err1)>ex,:) = [];
y_tmp1(abs(err1)>ex,:) = [];
p = polyfit(x_tmp1,y_tmp1,1);
ye = polyval(p,x_tmp1);
err = mean(y_tmp1-ye);
end