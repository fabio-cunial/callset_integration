FONT_SIZE=14;

supported=load('af_supported.txt');
unsupported=load('af_unsupported.txt');

% --------------------------------- N HAPS -------------------------------------
figure(1);
[y,x]=hist(supported(:,1),100); y=y./sum(y); 
subplot(1,2,1); hold on; plot(x,y,'.-');
subplot(1,2,2); hold on; y=cumsum(y); plot(x,y,'.-'); 
[y,x]=hist(unsupported(:,1),100); y=y./sum(y);
subplot(1,2,1); hold on; plot(x,y,'.-');
subplot(1,2,2); hold on; y=cumsum(y); plot(x,y,'.-'); 
subplot(1,2,1); set(gca,'fontsize',FONT_SIZE); xlabel('# haplotypes'); ylabel('% calls'); legend('supported','unsupported','fontsize',FONT_SIZE); axis square; grid on;
subplot(1,2,2); set(gca,'fontsize',FONT_SIZE); xlabel('# haplotypes'); ylabel('% calls with #haplotypes \leq'); axis square; grid on;

% ---------------------------------- AF ----------------------------------------
figure(2);
[y,x]=hist(supported(:,2),47); y=y./sum(y); 
subplot(1,2,1); hold on; plot(x,y,'.-');
subplot(1,2,2); hold on; y=cumsum(y); plot(x,y,'.-'); 
[y,x]=hist(unsupported(:,2),47); y=y./sum(y);
subplot(1,2,1); hold on; plot(x,y,'.-');
subplot(1,2,2); hold on; y=cumsum(y); plot(x,y,'.-'); 
subplot(1,2,1); set(gca,'fontsize',FONT_SIZE); xlabel('AF'); ylabel('% calls'); legend('supported','unsupported','fontsize',FONT_SIZE); axis square; grid on; axis([0,1,0,1]);
subplot(1,2,2); set(gca,'fontsize',FONT_SIZE); xlabel('AF'); ylabel('% calls with AF \leq'); axis square; grid on; axis([0,1,0,1]);

% ------------------------------ N SAMPLES -------------------------------------
figure(3);
[y,x]=hist(supported(:,3),47); y=y./sum(y); 
subplot(1,2,1); hold on; plot(x,y,'.-');
subplot(1,2,2); hold on; y=cumsum(y); plot(x,y,'.-'); 
[y,x]=hist(unsupported(:,3),47); y=y./sum(y);
subplot(1,2,1); hold on; plot(x,y,'.-');
subplot(1,2,2); hold on; y=cumsum(y); plot(x,y,'.-'); 
subplot(1,2,1); set(gca,'fontsize',FONT_SIZE); xlabel('# samples'); ylabel('% calls'); legend('supported','unsupported','fontsize',FONT_SIZE); axis square; grid on;
subplot(1,2,2); set(gca,'fontsize',FONT_SIZE); xlabel('# samples'); ylabel('% calls with #samples \leq'); axis square; grid on;
