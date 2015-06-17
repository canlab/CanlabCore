c = [1 -1 0]'
xcon = X * c;


for i = 1:1000
    e = mvnrnd(zeros(size(X,1),1),1);
    e = Vi * e;
    ppp(:,i) = psd(e);
end
pmean = mean(ppp,2);

figure; subplot(1,2,1);hold on; plot(psd(xcon),'b');
hold on;plot(pmean,'r');

subplot(1,2,2); h = cdfplot(psd(xcon)); hold on; h2 = cdfplot(psd(pmean)); set(h2,'Color','r')
legend({'Design' 'Error'})