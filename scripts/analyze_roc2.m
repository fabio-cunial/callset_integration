ROOT_DIR='/Users/fcunial/Downloads/sniffles-gt';
TR_STATUS=2;  % 0=no filter; 1=large overlap with the TR track; 2=small overlap with the TR track; 3=large overlap with the non-TR track.
TR_STATUS_LABELS={'All','TR large','TR small','Non-TR large'};


SVLEN_BINS={'100','500','2500','5000','10000','20000','50000'};
N_PANELS=7;
TAG_NAME='DV';
TAG_NAME_SNIFFLES='DV';
TAG_NAME_CUTESV='DV';
TAG_NAME_KANPIG='AD';
TAG_NAME_SVJEDIGRAPH='AD';

FONT_SIZE=18;
FONT_SIZE_LEGEND=12;


figure(1);


% ------------------------------- ROC CURVES -----------------------------------

% All calls
subplot(2,4,1); hold on;
SNIFFLES_FORCE=load(sprintf('%s/%d_sniffles_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_SNIFFLES));
CUTESV_FORCE=load(sprintf('%s/%d_cutesv_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_CUTESV));
KANPIG=load(sprintf('%s/%d_kanpig_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_KANPIG));
%SVJEDIGRAPH=load(sprintf('%s/%d_svjedigraph_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_SVJEDIGRAPH));
plot(SNIFFLES_FORCE(:,2),SNIFFLES_FORCE(:,3),'.-');
plot(CUTESV_FORCE(:,2),CUTESV_FORCE(:,3),'.-');
plot(KANPIG(:,2),KANPIG(:,3),'.-');
%plot(SVJEDIGRAPH(:,2),SVJEDIGRAPH(:,3),'.-');
line([0,1],[0,1],'color','black');
axis square; grid on; xlabel('False positive rate','fontsize',FONT_SIZE); 
ylabel('True positive rate','fontsize',FONT_SIZE);
title(sprintf('%s - %s',TR_STATUS_LABELS{TR_STATUS+1},TAG_NAME),'fontsize',FONT_SIZE);

% By SVLEN
for i=[2:N_PANELS] 
    subplot(2,4,i); hold on;
    LENGTH=SVLEN_BINS{i-1};
    SNIFFLES_FORCE=load(sprintf('%s/%d_sniffles_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_SNIFFLES));
    CUTESV_FORCE=load(sprintf('%s/%d_cutesv_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_CUTESV));
    KANPIG=load(sprintf('%s/%d_kanpig_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_KANPIG));
    %SVJEDIGRAPH=load(sprintf('%s/%d_svjedigraph_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_SVJEDIGRAPH));
    plot(SNIFFLES_FORCE(:,2),SNIFFLES_FORCE(:,3),'.-');
    plot(CUTESV_FORCE(:,2),CUTESV_FORCE(:,3),'.-');
    plot(KANPIG(:,2),KANPIG(:,3),'.-');
    %plot(SVJEDIGRAPH(:,2),SVJEDIGRAPH(:,3),'.-');
    line([0,1],[0,1],'color','black');
    axis square; grid on; xlabel('False positive rate','fontsize',FONT_SIZE); 
    ylabel('True positive rate','fontsize',FONT_SIZE);
    title(sprintf('%s - %s - %s',TR_STATUS_LABELS{TR_STATUS+1},TAG_NAME,LENGTH),'fontsize',FONT_SIZE);
endfor


% --------------------------- SUPPORTING CALLERS -------------------------------

A=load(sprintf('%s/%d_callers_merged.log',ROOT_DIR,TR_STATUS));
N_CALLS_ALL=A(:,2); N_DEL_ALL=A(:,3); N_INS_ALL=A(:,4); 
N_CALLS_2=A(:,5); N_DEL_2=A(:,6); N_INS_2=A(:,7);
N_CALLS_3=A(:,8); N_DEL_3=A(:,9); N_INS_3=A(:,10);
N_CALLS_PBSV_SNIFFLES=A(:,11); N_DEL_PBSV_SNIFFLES=A(:,12); N_INS_PBSV_SNIFFLES=A(:,13);
N_CALLS_PBSV_PAV=A(:,14); N_DEL_PBSV_PAV=A(:,15); N_INS_PBSV_PAV=A(:,16);
N_CALLS_SNIFFLES_PAV=A(:,17); N_DEL_SNIFFLES_PAV=A(:,18); N_INS_SNIFFLES_PAV=A(:,19);
N_CALLS_PBSV=A(:,20); N_DEL_PBSV=A(:,21); N_INS_PBSV=A(:,22); 
N_CALLS_SNIFFLES=A(:,23); N_DEL_SNIFFLES=A(:,24); N_INS_SNIFFLES=A(:,25); 
N_CALLS_PAV=A(:,26); N_DEL_PAV=A(:,27); N_INS_PAV=A(:,28); 

A=load(sprintf('%s/%d_callers_tp.log',ROOT_DIR,TR_STATUS));
N_TP=A(:,2); N_DEL_TP=A(:,3); N_INS_TP=A(:,4);
N_TP_2=A(:,5); N_DEL_TP_2=A(:,6); N_INS_TP_2=A(:,7);
N_TP_3=A(:,8); N_DEL_TP_3=A(:,9); N_INS_TP_3=A(:,10);
N_TP_PBSV_SNIFFLES=A(:,11); N_DEL_TP_PBSV_SNIFFLES=A(:,12); N_INS_TP_PBSV_SNIFFLES=A(:,13);
N_TP_PBSV_PAV=A(:,14); N_DEL_TP_PBSV_PAV=A(:,15); N_INS_TP_PBSV_PAV=A(:,16);
N_TP_SNIFFLES_PAV=A(:,17); N_DEL_TP_SNIFFLES_PAV=A(:,18); N_INS_TP_SNIFFLES_PAV=A(:,19);
N_TP_PBSV=A(:,20); N_DEL_TP_PBSV=A(:,21); N_INS_TP_PBSV=A(:,22); 
N_TP_SNIFFLES=A(:,23); N_DEL_TP_SNIFFLES=A(:,24); N_INS_TP_SNIFFLES=A(:,25); 
N_TP_PAV=A(:,26); N_DEL_TP_PAV=A(:,27); N_INS_TP_PAV=A(:,28); 

# >=2 CALLERS, all types.
TRUE_POSITIVES=N_TP_2;
FALSE_POSITIVES=N_CALLS_2-TRUE_POSITIVES;
TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
for i=[1:N_PANELS]
    subplot(2,4,i); hold on; plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'o');
endfor

# =3 CALLERS, all types.
TRUE_POSITIVES=N_TP_3;
FALSE_POSITIVES=N_CALLS_3-TRUE_POSITIVES;
TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
for i=[1:N_PANELS]
    subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'+');
endfor

# PBSV+SNIFFLES, all types.
TRUE_POSITIVES=N_TP_PBSV_SNIFFLES;
FALSE_POSITIVES=N_CALLS_PBSV_SNIFFLES-TRUE_POSITIVES;
TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
for i=[1:N_PANELS]
    subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'*');
endfor

# PBSV+PAV, all types.
TRUE_POSITIVES=N_TP_PBSV_PAV;
FALSE_POSITIVES=N_CALLS_PBSV_PAV-TRUE_POSITIVES;
TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
for i=[1:N_PANELS]
    subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'x');
endfor

# SNIFFLES+PAV, all types.
TRUE_POSITIVES=N_TP_SNIFFLES_PAV;
FALSE_POSITIVES=N_CALLS_SNIFFLES_PAV-TRUE_POSITIVES;
TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
for i=[1:N_PANELS]
    subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'s');
endfor

# PBSV, all types.
TRUE_POSITIVES=N_TP_PBSV;
FALSE_POSITIVES=N_CALLS_PBSV-TRUE_POSITIVES;
TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
for i=[1:N_PANELS]
    subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'^');
endfor

# SNIFFLES, all types.
TRUE_POSITIVES=N_TP_SNIFFLES;
FALSE_POSITIVES=N_CALLS_SNIFFLES-TRUE_POSITIVES;
TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
for i=[1:N_PANELS]
    subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'>');
endfor

# PAV, all types.
TRUE_POSITIVES=N_TP_PAV;
FALSE_POSITIVES=N_CALLS_PAV-TRUE_POSITIVES;
TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
for i=[1:N_PANELS]
    subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'v');
endfor
    

%legend('sniffles force','cutesv force','kanpig','svjedigraph','random','\geq2','=3','pbsv+snf','pbsv+pav','snf+pav','pbsv','sniffles','pav','location','southeast','fontsize',FONT_SIZE_LEGEND);
legend('sniffles force','cutesv force','kanpig','random','\geq2','=3','pbsv+snf','pbsv+pav','snf+pav','pbsv','sniffles','pav','location','southeast','fontsize',FONT_SIZE_LEGEND);
