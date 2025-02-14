% A1=load('/Users/fcunial/Downloads/TestILPCompression/latest/scip_easy_samples/runtimes-1.txt');
% A2=load('/Users/fcunial/Downloads/TestILPCompression/latest/scip_easy_samples/runtimes-2.txt');
% A3=load('/Users/fcunial/Downloads/TestILPCompression/latest/scip_easy_samples/runtimes-3.txt');
% A4=load('/Users/fcunial/Downloads/TestILPCompression/latest/scip_easy_samples/runtimes-4.txt');
% A=[A1;A2;A3;A4];

A=load('/Users/fcunial/Downloads/TestILPCompression/matrix.txt');

ASSUME_TIMEOUT_MS=15*60*1000;  % 15 mins


% WOULD_TIMEOUT_UNCOMPRESSED=length( find(A(:,7)>ASSUME_TIMEOUT_MS) )  % d
% WOULD_TIMEOUT_COMPRESSED=length( find(A(:,1)>ASSUME_TIMEOUT_MS) )  % d
%
% WOULD_TIMEOUT_UNCOMPRESSED=length( find(A(:,8)>ASSUME_TIMEOUT_MS) )  % n|d
% WOULD_TIMEOUT_COMPRESSED=length( find(A(:,3)>ASSUME_TIMEOUT_MS) )  % n|d

WOULD_TIMEOUT_UNCOMPRESSED=length( find(A(:,9)>ASSUME_TIMEOUT_MS) )  % d+n
WOULD_TIMEOUT_COMPRESSED=length( find(A(:,5)>ASSUME_TIMEOUT_MS) )  % d+n




figure(100);
B=[sum(A(:,7))./1000,sum(A(:,1))./1000,    sum(A(:,9))./1000,sum(A(:,5))./1000];
bar(B);
xticklabels({'d old','d new',    'd+n old','d+n new'});
ylabel('seconds'); title('Total runtime over all windows');
grid on;
set(gca,'fontsize',18);







figure(1);
B=[sum(A(:,7))./1000,sum(A(:,1))./1000,sum(A(:,2))./1000,    sum(A(:,8))./1000,sum(A(:,3))./1000,sum(A(:,4))./1000,    sum(A(:,9))./1000,sum(A(:,5))./1000,sum(A(:,6))./1000];
bar(B);
xticklabels({'d old','d new','d comp',    'n|d old','n|d new','n|d comp',    'd+n old','d+n new','d+n comp'});
ylabel('seconds'); title('Total runtime over all windows');
grid on;
set(gca,'fontsize',18);








figure(2);
B=[sum(A(:,13))./1000,sum(A(:,14))./1000,sum(A(:,15))./1000,sum(A(:,16))./1000,sum(A(:,17))./1000,sum(A(:,18))./1000,sum(A(:,19))./1000 ];
bar(B);
xticklabels({'mandatory','hapsGlobal','hapsLocal','reads','samples','hasLargeWeight','easySamples'});
ylabel('Total d+n runtime over all windows (seconds)'); title('Compression steps');
grid on;
set(gca,'fontsize',18);








figure(3);
nBins=100;

[y,x]=hist(A(:,1),nBins);
[yPrime,xPrime]=hist(A(:,7),nBins);
subplot(1,3,1); hold on;
loglog(x,y);
loglog(xPrime,yPrime);
title('d'); xlabel('milliseconds'); ylabel('n. windows');
axis square; grid on; set(gca,'fontsize',18);
legend('compressed','original');

[y,x]=hist(A(:,3),nBins);
[yPrime,xPrime]=hist(A(:,8),nBins);
subplot(1,3,2); hold on;
loglog(x,y);
loglog(xPrime,yPrime);
title('n given d'); xlabel('milliseconds'); ylabel('n. windows');
axis square; grid on; set(gca,'fontsize',18);

[y,x]=hist(A(:,5),nBins);
[yPrime,xPrime]=hist(A(:,9),nBins);
subplot(1,3,3); hold on;
loglog(x,y);
loglog(xPrime,yPrime);
title('d plus n'); xlabel('milliseconds'); ylabel('n. windows');
axis square; grid on; set(gca,'fontsize',18);


figure(4);

[y,x]=hist(A(:,7)./A(:,1),nBins);
subplot(1,3,1); hold on;
loglog(x,y);
title('d'); xlabel('speedup'); ylabel('n. windows');
axis square; grid on; set(gca,'fontsize',18);

[y,x]=hist(A(:,8)./A(:,3),nBins);
subplot(1,3,2); hold on;
loglog(x,y);
title('n given d'); xlabel('speedup'); ylabel('n. windows');
axis square; grid on; set(gca,'fontsize',18);

[y,x]=hist(A(:,9)./A(:,5),nBins);
subplot(1,3,3); hold on;
loglog(x,y);
title('d plus n'); xlabel('speedup'); ylabel('n. windows');
axis square; grid on; set(gca,'fontsize',18);


figure(5);

subplot(1,3,1);
plot(A(:,10),A(:,7)./A(:,1),'.');
title('d'); ylabel('speedup'); xlabel('n. egdes (original)');
axis square; grid on; set(gca,'fontsize',18);


subplot(1,3,2);
plot(A(:,11),A(:,8)./A(:,3),'.');
title('n given d'); ylabel('speedup'); xlabel('n. egdes (original)');
axis square; grid on; set(gca,'fontsize',18);


subplot(1,3,3);
plot(A(:,9),A(:,9)./A(:,5),'.');
title('d plus n'); ylabel('speedup'); xlabel('n. egdes (original)');
axis square; grid on; set(gca,'fontsize',18);
