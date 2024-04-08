FONT_SIZE=14;


% ----------------------------- SVLEN HISTOGRAM --------------------------------
lengths=[0:100:50000];
figure(1);

subplot(2,4,1); hold on;
A=load("./grch38/sniffles07/lengthHistogram.txt");
plot(lengths,A(:,1),'.-r'); plot(lengths,A(:,2),'.-b');
axis([0,6500,0,200000]); xlabel('SVLEN','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, sniffles 0.7','fontsize',FONT_SIZE);

subplot(2,4,2); hold on;
A=load("./grch38/sniffles09/lengthHistogram.txt");
plot(lengths,A(:,1),'.-r'); plot(lengths,A(:,2),'.-b');
axis([0,6500,0,200000]); xlabel('SVLEN','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, sniffles 0.9','fontsize',FONT_SIZE);

subplot(2,4,3); hold on;
A=load("./grch38/kanpig07/lengthHistogram.txt");
plot(lengths,A(:,1),'.-r'); plot(lengths,A(:,2),'.-b');
axis([0,6500,0,200000]); xlabel('SVLEN','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, kanpig 0.7','fontsize',FONT_SIZE);

subplot(2,4,4); hold on;
A=load("./grch38/kanpig09/lengthHistogram.txt");
plot(lengths,A(:,1),'.-r'); plot(lengths,A(:,2),'.-b');
axis([0,6500,0,200000]); xlabel('SVLEN','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, kanpig 0.9','fontsize',FONT_SIZE);

subplot(2,4,5); hold on;
A=load("./chm13/sniffles07/lengthHistogram.txt");
plot(lengths,A(:,1),'.-r'); plot(lengths,A(:,2),'.-b');
axis([0,6500,0,200000]); xlabel('SVLEN','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, sniffles 0.7','fontsize',FONT_SIZE);

subplot(2,4,6); hold on;
A=load("./chm13/sniffles09/lengthHistogram.txt");
plot(lengths,A(:,1),'.-r'); plot(lengths,A(:,2),'.-b');
axis([0,6500,0,200000]); xlabel('SVLEN','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, sniffles 0.9','fontsize',FONT_SIZE);

subplot(2,4,7); hold on;
A=load("./chm13/kanpig07/lengthHistogram.txt");
plot(lengths,A(:,1),'.-r'); plot(lengths,A(:,2),'.-b');
axis([0,6500,0,200000]); xlabel('SVLEN','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, kanpig 0.7','fontsize',FONT_SIZE);

subplot(2,4,8); hold on;
A=load("./chm13/kanpig09/lengthHistogram.txt");
plot(lengths,A(:,1),'.-r'); plot(lengths,A(:,2),'.-b');
axis([0,6500,0,200000]); xlabel('SVLEN','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, kanpig 0.9','fontsize',FONT_SIZE);




% ------------------------------- AF HISTOGRAM ---------------------------------
figure(2);

subplot(2,4,1); hold on;
A=load("./grch38/sniffles07/af.txt");
semilogy(A,'.-b');
xlabel('allele frequency','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, sniffles 0.7','fontsize',FONT_SIZE);

subplot(2,4,2); hold on;
A=load("./grch38/sniffles09/af.txt");
semilogy(A,'.-b');
xlabel('allele frequency','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, sniffles 0.9','fontsize',FONT_SIZE);

subplot(2,4,3); hold on;
A=load("./grch38/kanpig07/af.txt");
semilogy(A,'.-b');
xlabel('allele frequency','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, kanpig 0.7','fontsize',FONT_SIZE);

subplot(2,4,4); hold on;
A=load("./grch38/kanpig09/af.txt");
semilogy(A,'.-b');
xlabel('allele frequency','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, kanpig 0.9','fontsize',FONT_SIZE);

subplot(2,4,5); hold on;
A=load("./chm13/sniffles07/af.txt");
semilogy(A,'.-b');
xlabel('allele frequency','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, sniffles 0.7','fontsize',FONT_SIZE);

subplot(2,4,6); hold on;
A=load("./chm13/sniffles09/af.txt");
semilogy(A,'.-b');
xlabel('allele frequency','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, sniffles 0.9','fontsize',FONT_SIZE);

subplot(2,4,7); hold on;
A=load("./chm13/kanpig07/af.txt");
semilogy(A,'.-b');
xlabel('allele frequency','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, kanpig 0.7','fontsize',FONT_SIZE);

subplot(2,4,8); hold on;
A=load("./chm13/kanpig09/af.txt");
semilogy(A,'.-b');
xlabel('allele frequency','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, kanpig 0.9','fontsize',FONT_SIZE);




% ------------------------- INDIVIDUALS HISTOGRAM ------------------------------
DELTA=0.2;
figure(3);

subplot(2,4,1); hold on;
A=load("./grch38/sniffles07/individualHistogram.txt");
[nRows,nColumns]=size(A);
x=1-DELTA/2+rand(1,nRows)*DELTA; y=2-DELTA/2+rand(1,nRows)*DELTA;
plot(x,A(:,1),'or'); plot(y,A(:,2),'ob');
xlabel('SVTYPE','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, sniffles 0.7','fontsize',FONT_SIZE);

subplot(2,4,2); hold on;
A=load("./grch38/sniffles09/individualHistogram.txt");
[nRows,nColumns]=size(A);
plot(x,A(:,1),'or'); plot(y,A(:,2),'ob');
xlabel('SVTYPE','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, sniffles 0.9','fontsize',FONT_SIZE);

subplot(2,4,3); hold on;
A=load("./grch38/kanpig07/individualHistogram.txt");
[nRows,nColumns]=size(A);
plot(x,A(:,1),'or'); plot(y,A(:,2),'ob');
xlabel('SVTYPE','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, kanpig 0.7','fontsize',FONT_SIZE);

subplot(2,4,4); hold on;
A=load("./grch38/kanpig09/individualHistogram.txt");
[nRows,nColumns]=size(A);
plot(x,A(:,1),'or'); plot(y,A(:,2),'ob');
xlabel('SVTYPE','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('GRCh38, kanpig 0.9','fontsize',FONT_SIZE);

subplot(2,4,5); hold on;
A=load("./chm13/sniffles07/individualHistogram.txt");
[nRows,nColumns]=size(A);
x=1-DELTA/2+rand(1,nRows)*DELTA; y=2-DELTA/2+rand(1,nRows)*DELTA;
plot(x,A(:,1),'or'); plot(y,A(:,2),'ob');
xlabel('SVTYPE','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, sniffles 0.7','fontsize',FONT_SIZE);

subplot(2,4,6); hold on;
A=load("./chm13/sniffles09/individualHistogram.txt");
[nRows,nColumns]=size(A);
plot(x,A(:,1),'or'); plot(y,A(:,2),'ob');
xlabel('SVTYPE','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, sniffles 0.9','fontsize',FONT_SIZE);

subplot(2,4,7); hold on;
A=load("./chm13/kanpig07/individualHistogram.txt");
[nRows,nColumns]=size(A);
plot(x,A(:,1),'or'); plot(y,A(:,2),'ob');
xlabel('SVTYPE','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, kanpig 0.7','fontsize',FONT_SIZE);

subplot(2,4,8); hold on;
A=load("./chm13/kanpig09/individualHistogram.txt");
[nRows,nColumns]=size(A);
plot(x,A(:,1),'or'); plot(y,A(:,2),'ob');
xlabel('SVTYPE','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE); axis square; grid on;
title('CHM13, kanpig 0.9','fontsize',FONT_SIZE);
