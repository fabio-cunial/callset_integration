A=load('table_compressionTimes.txt');



figure(1);
B=[sum(A(:,7))./1000,sum(A(:,1))./1000,sum(A(:,2))./1000,    sum(A(:,8))./1000,sum(A(:,3))./1000,sum(A(:,4))./1000,    sum(A(:,9))./1000,sum(A(:,5))./1000,sum(A(:,6))./1000];
bar(B);
xticklabels({'d old','d new','d comp',    'n|d old','n|d new','n|d comp',    'd+n old','d+n new','d+n comp'});
ylabel('seconds'); title('Total runtime over all windows');
grid on;
set(gca,'fontsize',18);



figure(2);
B=[sum(A(:,13))./1000,sum(A(:,14))./1000,sum(A(:,15))./1000,sum(A(:,16))./1000,sum(A(:,17))./1000];
bar(B);
xticklabels({'mandatory','hapsGlobal','hapsLocal','reads','samples'});
ylabel('Total runtime over all windows (seconds)'); title('Compression steps');
grid on;
set(gca,'fontsize',18);








% figure(2);
% nBins=100;
%
% [y,x]=hist(A(:,1),nBins);
% [yPrime,xPrime]=hist(A(:,4),nBins);
% subplot(1,3,1); hold on;
% loglog(x,y);
% loglog(xPrime,yPrime);
% title('d'); xlabel('milliseconds'); ylabel('n. windows');
% axis square; grid on; set(gca,'fontsize',18);
% legend('compressed','original');
%
% [y,x]=hist(A(:,2),nBins);
% [yPrime,xPrime]=hist(A(:,5),nBins);
% subplot(1,3,2); hold on;
% loglog(x,y);
% loglog(xPrime,yPrime);
% title('n given d'); xlabel('milliseconds'); ylabel('n. windows');
% axis square; grid on; set(gca,'fontsize',18);
%
% [y,x]=hist(A(:,3),nBins);
% [yPrime,xPrime]=hist(A(:,6),nBins);
% subplot(1,3,3); hold on;
% loglog(x,y);
% loglog(xPrime,yPrime);
% title('d plus n'); xlabel('milliseconds'); ylabel('n. windows');
% axis square; grid on; set(gca,'fontsize',18);
%
%
% figure(3);
%
% [y,x]=hist(A(:,4)./A(:,1),nBins);
% subplot(1,3,1); hold on;
% loglog(x,y);
% title('d'); xlabel('speedup'); ylabel('n. windows');
% axis square; grid on; set(gca,'fontsize',18);
%
% [y,x]=hist(A(:,5)./A(:,2),nBins);
% subplot(1,3,2); hold on;
% loglog(x,y);
% title('n given d'); xlabel('speedup'); ylabel('n. windows');
% axis square; grid on; set(gca,'fontsize',18);
%
% [y,x]=hist(A(:,6)./A(:,3),nBins);
% subplot(1,3,3); hold on;
% loglog(x,y);
% title('d plus n'); xlabel('speedup'); ylabel('n. windows');
% axis square; grid on; set(gca,'fontsize',18);
%
%
% figure(4);
%
% subplot(1,3,1);
% plot(A(:,6),A(:,4)./A(:,1),'.');
% title('d'); ylabel('speedup'); xlabel('n. egdes (original)');
% axis square; grid on; set(gca,'fontsize',18);
%
%
% subplot(1,3,2);
% plot(A(:,7),A(:,5)./A(:,2),'.');
% title('n given d'); ylabel('speedup'); xlabel('n. egdes (original)');
% axis square; grid on; set(gca,'fontsize',18);
%
%
% subplot(1,3,3);
% plot(A(:,8),A(:,6)./A(:,3),'.');
% title('d plus n'); ylabel('speedup'); xlabel('n. egdes (original)');
% axis square; grid on; set(gca,'fontsize',18);
