# Analyzes the output of $VCFstats.java$.
#
FONT_SIZE=14;
COLOR_UNFILTERED=[147,232,249]./256;
COLOR_09=[17,205,242]./256;
COLOR_07=[9,145,172]./256;
N_SAMPLES=47;

TITLE='GRCh38, sniffles';


figure(1);
subplot(1,3,1); hold on;
A=load('./hprc/grch38/support_unfiltered.txt');
h=plot3(A(:,1),A(:,2),A(:,3),'.'); set(h,'color',COLOR_UNFILTERED);
axis square; grid on; title(TITLE,'fontsize',FONT_SIZE);
xlabel('pbsv','fontsize',FONT_SIZE); ylabel('sniffles','fontsize',FONT_SIZE); zlabel('pav','fontsize',FONT_SIZE);

figure(2);
subplot(1,3,1); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,1)+1,A(i,2)+1)=B(A(i,1)+1,A(i,2)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('No filtering','fontsize',FONT_SIZE); ylabel('pbsv','fontsize',FONT_SIZE); xlabel('sniffles','fontsize',FONT_SIZE);

figure(3);
subplot(1,3,1); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,1)+1,A(i,3)+1)=B(A(i,1)+1,A(i,3)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('No filtering','fontsize',FONT_SIZE); ylabel('pbsv','fontsize',FONT_SIZE); xlabel('pav','fontsize',FONT_SIZE);

figure(4);
subplot(1,3,1); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,2)+1,A(i,3)+1)=B(A(i,2)+1,A(i,3)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('No filtering','fontsize',FONT_SIZE); ylabel('sniffles','fontsize',FONT_SIZE); xlabel('pav','fontsize',FONT_SIZE);


figure(1);
A=load('./hprc/grch38/af_sniffles_unfiltered.txt');
subplot(1,3,2); hold on;
h=plot([0:length(A)-1],A,'.-'); set(h,'color',COLOR_UNFILTERED);
xlabel('AF','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on;

A=load('./hprc/grch38/sf_sniffles_unfiltered.txt');
subplot(1,3,3); hold on;
h=plot([0:length(A)-1],A,'.-'); set(h,'color',COLOR_UNFILTERED);
xlabel('SF','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on;



figure(1);
subplot(1,3,1); hold on;
A=load('./hprc/grch38/support_sniffles_09.txt');
h=plot3(A(:,1),A(:,2),A(:,3),'.'); set(h,'color',COLOR_09);
axis square; grid on; title(TITLE,'fontsize',FONT_SIZE);
xlabel('pbsv','fontsize',FONT_SIZE); ylabel('sniffles','fontsize',FONT_SIZE); zlabel('pav','fontsize',FONT_SIZE);

figure(2);
subplot(1,3,2); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,1)+1,A(i,2)+1)=B(A(i,1)+1,A(i,2)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('TPR=0.9','fontsize',FONT_SIZE); ylabel('pbsv','fontsize',FONT_SIZE); xlabel('sniffles','fontsize',FONT_SIZE);

figure(3);
subplot(1,3,2); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,1)+1,A(i,3)+1)=B(A(i,1)+1,A(i,3)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('TPR=0.9','fontsize',FONT_SIZE); ylabel('pbsv','fontsize',FONT_SIZE); xlabel('pav','fontsize',FONT_SIZE);

figure(4);
subplot(1,3,2); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,2)+1,A(i,3)+1)=B(A(i,2)+1,A(i,3)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('TPR=0.9','fontsize',FONT_SIZE); ylabel('sniffles','fontsize',FONT_SIZE); xlabel('pav','fontsize',FONT_SIZE);


figure(1);
A=load('./hprc/grch38/af_sniffles_09.txt');
subplot(1,3,2); hold on;
h=plot([0:length(A)-1],A,'.-'); set(h,'color',COLOR_09);
xlabel('AF','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on;

A=load('./hprc/grch38/sf_sniffles_09.txt');
subplot(1,3,3); hold on;
h=plot([0:length(A)-1],A,'.-'); set(h,'color',COLOR_09);
xlabel('SF','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on;



figure(1);
subplot(1,3,1); hold on;
A=load('./hprc/grch38/support_sniffles_07.txt');
h=plot3(A(:,1),A(:,2),A(:,3),'.'); set(h,'color',COLOR_07);
axis square; grid on; title(TITLE,'fontsize',FONT_SIZE);
xlabel('pbsv','fontsize',FONT_SIZE); ylabel('sniffles','fontsize',FONT_SIZE); zlabel('pav','fontsize',FONT_SIZE);

figure(2)
subplot(1,3,3); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,1)+1,A(i,2)+1)=B(A(i,1)+1,A(i,2)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('TPR=0.7','fontsize',FONT_SIZE); ylabel('pbsv','fontsize',FONT_SIZE); xlabel('sniffles','fontsize',FONT_SIZE);

figure(3)
subplot(1,3,3); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,1)+1,A(i,3)+1)=B(A(i,1)+1,A(i,3)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('TPR=0.7','fontsize',FONT_SIZE); ylabel('pbsv','fontsize',FONT_SIZE); xlabel('pav','fontsize',FONT_SIZE);

figure(4)
subplot(1,3,3); hold on;
B=zeros(N_SAMPLES+1,N_SAMPLES+1);
[nRows,nColums]=size(A);
for i=[1:nRows]
    B(A(i,2)+1,A(i,3)+1)=B(A(i,2)+1,A(i,3)+1)+1;
endfor
colormap('hot'); imagesc(log10(B)); caxis([0,5]);
axis square; grid on; title('TPR=0.7','fontsize',FONT_SIZE); ylabel('sniffles','fontsize',FONT_SIZE); xlabel('pav','fontsize',FONT_SIZE);


figure(1);
A=load('./hprc/grch38/af_sniffles_07.txt');
subplot(1,3,2); hold on;
h=plot([0:length(A)-1],A,'.-'); set(h,'color',COLOR_07);
xlabel('AF','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on;

A=load('./hprc/grch38/sf_sniffles_07.txt');
subplot(1,3,3); hold on;
h=plot([0:length(A)-1],A,'.-'); set(h,'color',COLOR_07);
xlabel('SF','fontsize',FONT_SIZE); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on;


figure(5);
N_CALLERS=3;
A1=load('./hprc/grch38/rareCalls_sniffles_unfiltered.txt');
A2=load('./hprc/grch38/rareCalls_sniffles_09.txt');
A3=load('./hprc/grch38/rareCalls_sniffles_07.txt');

subplot(1,4,1); hold on;
B=[  A1(1,2:N_CALLERS+1),A1(1,N_CALLERS+2)-sum(A1(1,2:N_CALLERS+1));  A2(1,2:N_CALLERS+1),A2(1,N_CALLERS+2)-sum(A2(1,2:N_CALLERS+1));  A3(1,2:N_CALLERS+1),A3(1,N_CALLERS+2)-sum(A3(1,2:N_CALLERS+1))  ];
bar(B,'stacked'); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on; 
legend('pbsv only','sniffles only','pav only','rest');
title('Calls that occur in no sample','fontsize',FONT_SIZE);
xticks([1,2,3]); xticklabels({'unfiltered','0.9','0.7'});

subplot(1,4,2); hold on;
B=[  A1(2,2:N_CALLERS+1),A1(2,N_CALLERS+2)-sum(A1(2,2:N_CALLERS+1));  A2(2,2:N_CALLERS+1),A2(2,N_CALLERS+2)-sum(A2(2,2:N_CALLERS+1));  A3(2,2:N_CALLERS+1),A3(2,N_CALLERS+2)-sum(A3(2,2:N_CALLERS+1))  ];
bar(B,'stacked'); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on; 
%legend('pbsv only','sniffles only','pav only','rest');
title('Calls that occur in one sample','fontsize',FONT_SIZE);
xticks([1,2,3]); xticklabels({'unfiltered','0.9','0.7'});

subplot(1,4,3); hold on;
B=[  A1(3,2:N_CALLERS+1),A1(3,N_CALLERS+2)-sum(A1(3,2:N_CALLERS+1));  A2(3,2:N_CALLERS+1),A2(3,N_CALLERS+2)-sum(A2(3,2:N_CALLERS+1));  A3(3,2:N_CALLERS+1),A3(3,N_CALLERS+2)-sum(A3(3,2:N_CALLERS+1))  ];
bar(B,'stacked'); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on; 
%legend('pbsv only','sniffles only','pav only','rest');
title('Calls that occur in two samples','fontsize',FONT_SIZE);
xticks([1,2,3]); xticklabels({'unfiltered','0.9','0.7'});

subplot(1,4,4); hold on;
B=[  A1(4,2:N_CALLERS+1),A1(4,N_CALLERS+2)-sum(A1(4,2:N_CALLERS+1));  A2(4,2:N_CALLERS+1),A2(4,N_CALLERS+2)-sum(A2(4,2:N_CALLERS+1));  A3(4,2:N_CALLERS+1),A3(4,N_CALLERS+2)-sum(A3(4,2:N_CALLERS+1))  ];
bar(B,'stacked'); ylabel('# records','fontsize',FONT_SIZE);
axis square; grid on; 
%legend('pbsv only','sniffles only','pav only','rest');
title('Calls that occur in three samples','fontsize',FONT_SIZE);
xticks([1,2,3]); xticklabels({'unfiltered','0.9','0.7'});