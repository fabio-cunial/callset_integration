ROOT_DIR='/Users/fcunial/Downloads/sniffles-gt/HPRC';
TR_STATUS=0;  % 0=no filter; 1=large overlap with the TR track; 2=small overlap with the TR track; 3=large overlap with the non-TR track.
TR_STATUS_LABELS={'All','TR 0.9','TR 0.1','Non-TR 0.9'};


SVLEN_BINS={'100','500','2500','5000','10000','20000','50000'};
SVLEN_BIN_LABELS={'[50..100)','[100..500)','[500..2500)','[2500..5k)','[5k..10k)','[10k..20k)','[20k..50k)'};
SAMPLES={'HG002','HG005','HG00621','HG00733','HG00735','HG00741','HG01071','HG01106','HG01109','HG01123','HG01243','HG01358','HG01891','HG01928','HG01952','HG01978','HG02055','HG02080','HG02109','HG02145','HG02148','HG02486','HG02559','HG02572','HG02622','HG02630','HG02717','HG02818','HG02886','HG03098','HG03453','HG03486','HG03492','HG03540','HG03579','NA18906','NA19240','NA21309'};
N_PANELS=7;
LINE_WIDTH_SAMS=2;

COLOR_SNIFFLES=[0,114,189]./256;
COLOR_CUTESV=[217,83,25]./256;
COLOR_KANPIG=[237,177,32]./256;
COLOR_GEQ2=[119,172,48]./256;
COLOR_THREE=[77,190,238]./256;
COLOR_PBSV_SNIFFLES=[163,24,50]./256;
COLOR_PBSV_PAV=[0,114,189]./256;
COLOR_SNIFFLES_PAV=[217,83,25]./256;
COLOR_PBSV=[237,177,32]./256;
COLOR_SNIFFLES_ONLY=[126,47,142]./256;
COLOR_PAV=[119,172,48]./256;

% TAG_NAME='DV';
% TAG_NAME_SNIFFLES='DV';
% TAG_NAME_CUTESV='DV';
% TAG_NAME_KANPIG='AD';

TAG_NAME='GQ';
TAG_NAME_SNIFFLES='GQ';
TAG_NAME_CUTESV='GQ';
TAG_NAME_KANPIG='SQ';

FONT_SIZE=18;
FONT_SIZE_LEGEND=12;





% ------------------------------- MAIN SLIDE -----------------------------------
TR_STATUS_PRIME=TR_STATUS;  % Backup copy to be restored later
figure(2);

% ---------------------------- ROC CURVES NON-TR -------------------------------
TR_STATUS=3; subplot(1,2,1); hold on;

for i=[1:length(SAMPLES)]
    SNIFFLES_FORCE=load(sprintf('%s/%s_logs/%d_sniffles_roc_calls_%s_geq.log',ROOT_DIR,SAMPLES{i},TR_STATUS,TAG_NAME_SNIFFLES));
    CUTESV_FORCE=load(sprintf('%s/%s_logs/%d_cutesv_roc_calls_%s_geq.log',ROOT_DIR,SAMPLES{i},TR_STATUS,TAG_NAME_CUTESV));
    KANPIG=load(sprintf('%s/%s_logs/%d_kanpig_roc_calls_%s_geq.log',ROOT_DIR,SAMPLES{i},TR_STATUS,TAG_NAME_KANPIG));
    h=plot(SNIFFLES_FORCE(:,2),SNIFFLES_FORCE(:,3),'-'); set(h,'color',COLOR_SNIFFLES);
    h=plot(CUTESV_FORCE(:,2),CUTESV_FORCE(:,3),'-'); set(h,'color',COLOR_CUTESV);
    h=plot(KANPIG(:,2),KANPIG(:,3),'-'); set(h,'color',COLOR_KANPIG);
endfor
line([0,1],[0,1],'color','black');
axis square; grid on; xlabel('False positive rate','fontsize',FONT_SIZE); 
ylabel('True positive rate','fontsize',FONT_SIZE);
title(sprintf('%s, %s',TR_STATUS_LABELS{TR_STATUS+1},TAG_NAME),'fontsize',FONT_SIZE);

% ------------------------ SUPPORTING CALLERS NON-TR ---------------------------

for j=[1:length(SAMPLES)]
    A=load(sprintf('%s/%s_logs/%d_callers_merged.log',ROOT_DIR,SAMPLES{j},TR_STATUS));
    N_CALLS_ALL=A(:,2); N_DEL_ALL=A(:,3); N_INS_ALL=A(:,4); 
    N_CALLS_2=A(:,5); N_DEL_2=A(:,6); N_INS_2=A(:,7);
    N_CALLS_3=A(:,8); N_DEL_3=A(:,9); N_INS_3=A(:,10);
    N_CALLS_PBSV_SNIFFLES=A(:,11); N_DEL_PBSV_SNIFFLES=A(:,12); N_INS_PBSV_SNIFFLES=A(:,13);
    N_CALLS_PBSV_PAV=A(:,14); N_DEL_PBSV_PAV=A(:,15); N_INS_PBSV_PAV=A(:,16);
    N_CALLS_SNIFFLES_PAV=A(:,17); N_DEL_SNIFFLES_PAV=A(:,18); N_INS_SNIFFLES_PAV=A(:,19);
    N_CALLS_PBSV=A(:,20); N_DEL_PBSV=A(:,21); N_INS_PBSV=A(:,22); 
    N_CALLS_SNIFFLES=A(:,23); N_DEL_SNIFFLES=A(:,24); N_INS_SNIFFLES=A(:,25); 
    N_CALLS_PAV=A(:,26); N_DEL_PAV=A(:,27); N_INS_PAV=A(:,28); 

    A=load(sprintf('%s/%s_logs/%d_callers_tp.log',ROOT_DIR,SAMPLES{j},TR_STATUS));
    N_TP=A(:,2); N_DEL_TP=A(:,3); N_INS_TP=A(:,4);
    N_TP_2=A(:,5); N_DEL_TP_2=A(:,6); N_INS_TP_2=A(:,7);
    N_TP_3=A(:,8); N_DEL_TP_3=A(:,9); N_INS_TP_3=A(:,10);
    N_TP_PBSV_SNIFFLES=A(:,11); N_DEL_TP_PBSV_SNIFFLES=A(:,12); N_INS_TP_PBSV_SNIFFLES=A(:,13);
    N_TP_PBSV_PAV=A(:,14); N_DEL_TP_PBSV_PAV=A(:,15); N_INS_TP_PBSV_PAV=A(:,16);
    N_TP_SNIFFLES_PAV=A(:,17); N_DEL_TP_SNIFFLES_PAV=A(:,18); N_INS_TP_SNIFFLES_PAV=A(:,19);
    N_TP_PBSV=A(:,20); N_DEL_TP_PBSV=A(:,21); N_INS_TP_PBSV=A(:,22); 
    N_TP_SNIFFLES=A(:,23); N_DEL_TP_SNIFFLES=A(:,24); N_INS_TP_SNIFFLES=A(:,25); 
    N_TP_PAV=A(:,26); N_DEL_TP_PAV=A(:,27); N_INS_TP_PAV=A(:,28); 

    i=1;

    # >=2 CALLERS, all types.
    TRUE_POSITIVES=N_TP_2;
    FALSE_POSITIVES=N_CALLS_2-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'o'); set(h,'color',COLOR_GEQ2);

    # =3 CALLERS, all types.
    TRUE_POSITIVES=N_TP_3;
    FALSE_POSITIVES=N_CALLS_3-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'+'); set(h,'color',COLOR_THREE);

    # PBSV+SNIFFLES, all types.
    TRUE_POSITIVES=N_TP_PBSV_SNIFFLES;
    FALSE_POSITIVES=N_CALLS_PBSV_SNIFFLES-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'*'); set(h,'color',COLOR_PBSV_SNIFFLES);

    # PBSV+PAV, all types.
    TRUE_POSITIVES=N_TP_PBSV_PAV;
    FALSE_POSITIVES=N_CALLS_PBSV_PAV-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'x'); set(h,'color',COLOR_PBSV_PAV);

    # SNIFFLES+PAV, all types.
    TRUE_POSITIVES=N_TP_SNIFFLES_PAV;
    FALSE_POSITIVES=N_CALLS_SNIFFLES_PAV-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'s'); set(h,'color',COLOR_SNIFFLES_PAV);

    # PBSV, all types.
    TRUE_POSITIVES=N_TP_PBSV;
    FALSE_POSITIVES=N_CALLS_PBSV-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'^'); set(h,'color',COLOR_PBSV);

    # SNIFFLES, all types.
    TRUE_POSITIVES=N_TP_SNIFFLES;
    FALSE_POSITIVES=N_CALLS_SNIFFLES-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'>'); set(h,'color',COLOR_SNIFFLES_ONLY);

    # PAV, all types.
    TRUE_POSITIVES=N_TP_PAV;
    FALSE_POSITIVES=N_CALLS_PAV-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'v'); set(h,'color',COLOR_PAV);
endfor


% ------------------------------ ROC CURVES TR ---------------------------------

TR_STATUS=1; subplot(1,2,2); hold on;

for i=[1:length(SAMPLES)]
    SNIFFLES_FORCE=load(sprintf('%s/%s_logs/%d_sniffles_roc_calls_%s_geq.log',ROOT_DIR,SAMPLES{i},TR_STATUS,TAG_NAME_SNIFFLES));
    CUTESV_FORCE=load(sprintf('%s/%s_logs/%d_cutesv_roc_calls_%s_geq.log',ROOT_DIR,SAMPLES{i},TR_STATUS,TAG_NAME_CUTESV));
    KANPIG=load(sprintf('%s/%s_logs/%d_kanpig_roc_calls_%s_geq.log',ROOT_DIR,SAMPLES{i},TR_STATUS,TAG_NAME_KANPIG));
    h=plot(SNIFFLES_FORCE(:,2),SNIFFLES_FORCE(:,3),'-'); set(h,'color',COLOR_SNIFFLES);
    h=plot(CUTESV_FORCE(:,2),CUTESV_FORCE(:,3),'-'); set(h,'color',COLOR_CUTESV);
    h=plot(KANPIG(:,2),KANPIG(:,3),'-'); set(h,'color',COLOR_KANPIG);
endfor
line([0,1],[0,1],'color','black');
axis square; grid on; xlabel('False positive rate','fontsize',FONT_SIZE); 
ylabel('True positive rate','fontsize',FONT_SIZE);
title(sprintf('%s, %s',TR_STATUS_LABELS{TR_STATUS+1},TAG_NAME),'fontsize',FONT_SIZE);

% ------------------------ SUPPORTING CALLERS NON-TR ---------------------------

for j=[1:length(SAMPLES)]
    A=load(sprintf('%s/%s_logs/%d_callers_merged.log',ROOT_DIR,SAMPLES{j},TR_STATUS));
    N_CALLS_ALL=A(:,2); N_DEL_ALL=A(:,3); N_INS_ALL=A(:,4); 
    N_CALLS_2=A(:,5); N_DEL_2=A(:,6); N_INS_2=A(:,7);
    N_CALLS_3=A(:,8); N_DEL_3=A(:,9); N_INS_3=A(:,10);
    N_CALLS_PBSV_SNIFFLES=A(:,11); N_DEL_PBSV_SNIFFLES=A(:,12); N_INS_PBSV_SNIFFLES=A(:,13);
    N_CALLS_PBSV_PAV=A(:,14); N_DEL_PBSV_PAV=A(:,15); N_INS_PBSV_PAV=A(:,16);
    N_CALLS_SNIFFLES_PAV=A(:,17); N_DEL_SNIFFLES_PAV=A(:,18); N_INS_SNIFFLES_PAV=A(:,19);
    N_CALLS_PBSV=A(:,20); N_DEL_PBSV=A(:,21); N_INS_PBSV=A(:,22); 
    N_CALLS_SNIFFLES=A(:,23); N_DEL_SNIFFLES=A(:,24); N_INS_SNIFFLES=A(:,25); 
    N_CALLS_PAV=A(:,26); N_DEL_PAV=A(:,27); N_INS_PAV=A(:,28); 

    A=load(sprintf('%s/%s_logs/%d_callers_tp.log',ROOT_DIR,SAMPLES{j},TR_STATUS));
    N_TP=A(:,2); N_DEL_TP=A(:,3); N_INS_TP=A(:,4);
    N_TP_2=A(:,5); N_DEL_TP_2=A(:,6); N_INS_TP_2=A(:,7);
    N_TP_3=A(:,8); N_DEL_TP_3=A(:,9); N_INS_TP_3=A(:,10);
    N_TP_PBSV_SNIFFLES=A(:,11); N_DEL_TP_PBSV_SNIFFLES=A(:,12); N_INS_TP_PBSV_SNIFFLES=A(:,13);
    N_TP_PBSV_PAV=A(:,14); N_DEL_TP_PBSV_PAV=A(:,15); N_INS_TP_PBSV_PAV=A(:,16);
    N_TP_SNIFFLES_PAV=A(:,17); N_DEL_TP_SNIFFLES_PAV=A(:,18); N_INS_TP_SNIFFLES_PAV=A(:,19);
    N_TP_PBSV=A(:,20); N_DEL_TP_PBSV=A(:,21); N_INS_TP_PBSV=A(:,22); 
    N_TP_SNIFFLES=A(:,23); N_DEL_TP_SNIFFLES=A(:,24); N_INS_TP_SNIFFLES=A(:,25); 
    N_TP_PAV=A(:,26); N_DEL_TP_PAV=A(:,27); N_INS_TP_PAV=A(:,28); 

    i=1;

    # >=2 CALLERS, all types.
    TRUE_POSITIVES=N_TP_2;
    FALSE_POSITIVES=N_CALLS_2-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'o'); set(h,'color',COLOR_GEQ2);

    # =3 CALLERS, all types.
    TRUE_POSITIVES=N_TP_3;
    FALSE_POSITIVES=N_CALLS_3-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'+'); set(h,'color',COLOR_THREE);

    # PBSV+SNIFFLES, all types.
    TRUE_POSITIVES=N_TP_PBSV_SNIFFLES;
    FALSE_POSITIVES=N_CALLS_PBSV_SNIFFLES-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'*'); set(h,'color',COLOR_PBSV_SNIFFLES);

    # PBSV+PAV, all types.
    TRUE_POSITIVES=N_TP_PBSV_PAV;
    FALSE_POSITIVES=N_CALLS_PBSV_PAV-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'x'); set(h,'color',COLOR_PBSV_PAV);

    # SNIFFLES+PAV, all types.
    TRUE_POSITIVES=N_TP_SNIFFLES_PAV;
    FALSE_POSITIVES=N_CALLS_SNIFFLES_PAV-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'s'); set(h,'color',COLOR_SNIFFLES_PAV);

    # PBSV, all types.
    TRUE_POSITIVES=N_TP_PBSV;
    FALSE_POSITIVES=N_CALLS_PBSV-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'^'); set(h,'color',COLOR_PBSV);

    # SNIFFLES, all types.
    TRUE_POSITIVES=N_TP_SNIFFLES;
    FALSE_POSITIVES=N_CALLS_SNIFFLES-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'>'); set(h,'color',COLOR_SNIFFLES_ONLY);

    # PAV, all types.
    TRUE_POSITIVES=N_TP_PAV;
    FALSE_POSITIVES=N_CALLS_PAV-TRUE_POSITIVES;
    TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
    FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
    h=plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'v'); set(h,'color',COLOR_PAV);
endfor



% % -------------------------------- BY SVLEN ------------------------------------
% TR_STATUS=TR_STATUS_PRIME;  % Restoring from backup
%
% % ------------------------------- ROC CURVES -----------------------------------
% figure(1);
%
% % All calls
% subplot(2,4,1); hold on;
% SNIFFLES_FORCE=load(sprintf('%s/%d_sniffles_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_SNIFFLES));
% CUTESV_FORCE=load(sprintf('%s/%d_cutesv_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_CUTESV));
% KANPIG=load(sprintf('%s/%d_kanpig_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_KANPIG));
% %SVJEDIGRAPH=load(sprintf('%s/%d_svjedigraph_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_SVJEDIGRAPH));
% SAMS=load(sprintf('%s/%d_sams_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,TAG_NAME_SAMS));
% plot(SNIFFLES_FORCE(:,2),SNIFFLES_FORCE(:,3),'.-');
% plot(CUTESV_FORCE(:,2),CUTESV_FORCE(:,3),'.-');
% plot(KANPIG(:,2),KANPIG(:,3),'.-');
% %plot(SVJEDIGRAPH(:,2),SVJEDIGRAPH(:,3),'.-');
% h=plot(SAMS(:,2),SAMS(:,3),'-'); set(h,'linewidth',LINE_WIDTH_SAMS);
% line([0,1],[0,1],'color','black');
% axis square; grid on; xlabel('False positive rate','fontsize',FONT_SIZE);
% ylabel('True positive rate','fontsize',FONT_SIZE);
% title(sprintf('%s, %s',TR_STATUS_LABELS{TR_STATUS+1},TAG_NAME),'fontsize',FONT_SIZE);
%
% % By SVLEN
% for i=[2:N_PANELS]
%     subplot(2,4,i); hold on;
%     LENGTH=SVLEN_BINS{i-1};
%     SNIFFLES_FORCE=load(sprintf('%s/%d_sniffles_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_SNIFFLES));
%     CUTESV_FORCE=load(sprintf('%s/%d_cutesv_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_CUTESV));
%     KANPIG=load(sprintf('%s/%d_kanpig_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_KANPIG));
%     %SVJEDIGRAPH=load(sprintf('%s/%d_svjedigraph_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_SVJEDIGRAPH));
%     SAMS=load(sprintf('%s/%d_sams_%s_roc_calls_%s_geq.log',ROOT_DIR,TR_STATUS,LENGTH,TAG_NAME_SAMS));
%     plot(SNIFFLES_FORCE(:,2),SNIFFLES_FORCE(:,3),'.-');
%     plot(CUTESV_FORCE(:,2),CUTESV_FORCE(:,3),'.-');
%     plot(KANPIG(:,2),KANPIG(:,3),'.-');
%     %plot(SVJEDIGRAPH(:,2),SVJEDIGRAPH(:,3),'.-');
%     h=plot(SAMS(:,2),SAMS(:,3),'-'); set(h,'linewidth',LINE_WIDTH_SAMS);
%     line([0,1],[0,1],'color','black');
%     axis square; grid on; xlabel('False positive rate','fontsize',FONT_SIZE);
%     ylabel('True positive rate','fontsize',FONT_SIZE);
%     title(sprintf('%s, %s, %s',TR_STATUS_LABELS{TR_STATUS+1},SVLEN_BIN_LABELS{i-1},TAG_NAME),'fontsize',FONT_SIZE);
% endfor
%
%
% % --------------------------- SUPPORTING CALLERS -------------------------------
%
% A=load(sprintf('%s/%d_callers_merged.log',ROOT_DIR,TR_STATUS));
% N_CALLS_ALL=A(:,2); N_DEL_ALL=A(:,3); N_INS_ALL=A(:,4);
% N_CALLS_2=A(:,5); N_DEL_2=A(:,6); N_INS_2=A(:,7);
% N_CALLS_3=A(:,8); N_DEL_3=A(:,9); N_INS_3=A(:,10);
% N_CALLS_PBSV_SNIFFLES=A(:,11); N_DEL_PBSV_SNIFFLES=A(:,12); N_INS_PBSV_SNIFFLES=A(:,13);
% N_CALLS_PBSV_PAV=A(:,14); N_DEL_PBSV_PAV=A(:,15); N_INS_PBSV_PAV=A(:,16);
% N_CALLS_SNIFFLES_PAV=A(:,17); N_DEL_SNIFFLES_PAV=A(:,18); N_INS_SNIFFLES_PAV=A(:,19);
% N_CALLS_PBSV=A(:,20); N_DEL_PBSV=A(:,21); N_INS_PBSV=A(:,22);
% N_CALLS_SNIFFLES=A(:,23); N_DEL_SNIFFLES=A(:,24); N_INS_SNIFFLES=A(:,25);
% N_CALLS_PAV=A(:,26); N_DEL_PAV=A(:,27); N_INS_PAV=A(:,28);
%
% A=load(sprintf('%s/%d_callers_tp.log',ROOT_DIR,TR_STATUS));
% N_TP=A(:,2); N_DEL_TP=A(:,3); N_INS_TP=A(:,4);
% N_TP_2=A(:,5); N_DEL_TP_2=A(:,6); N_INS_TP_2=A(:,7);
% N_TP_3=A(:,8); N_DEL_TP_3=A(:,9); N_INS_TP_3=A(:,10);
% N_TP_PBSV_SNIFFLES=A(:,11); N_DEL_TP_PBSV_SNIFFLES=A(:,12); N_INS_TP_PBSV_SNIFFLES=A(:,13);
% N_TP_PBSV_PAV=A(:,14); N_DEL_TP_PBSV_PAV=A(:,15); N_INS_TP_PBSV_PAV=A(:,16);
% N_TP_SNIFFLES_PAV=A(:,17); N_DEL_TP_SNIFFLES_PAV=A(:,18); N_INS_TP_SNIFFLES_PAV=A(:,19);
% N_TP_PBSV=A(:,20); N_DEL_TP_PBSV=A(:,21); N_INS_TP_PBSV=A(:,22);
% N_TP_SNIFFLES=A(:,23); N_DEL_TP_SNIFFLES=A(:,24); N_INS_TP_SNIFFLES=A(:,25);
% N_TP_PAV=A(:,26); N_DEL_TP_PAV=A(:,27); N_INS_TP_PAV=A(:,28);
%
% # >=2 CALLERS, all types.
% TRUE_POSITIVES=N_TP_2;
% FALSE_POSITIVES=N_CALLS_2-TRUE_POSITIVES;
% TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
% FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
% for i=[1:N_PANELS]
%     subplot(2,4,i); hold on; plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'o');
% endfor
%
% # =3 CALLERS, all types.
% TRUE_POSITIVES=N_TP_3;
% FALSE_POSITIVES=N_CALLS_3-TRUE_POSITIVES;
% TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
% FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
% for i=[1:N_PANELS]
%     subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'+');
% endfor
%
% # PBSV+SNIFFLES, all types.
% TRUE_POSITIVES=N_TP_PBSV_SNIFFLES;
% FALSE_POSITIVES=N_CALLS_PBSV_SNIFFLES-TRUE_POSITIVES;
% TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
% FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
% for i=[1:N_PANELS]
%     subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'*');
% endfor
%
% # PBSV+PAV, all types.
% TRUE_POSITIVES=N_TP_PBSV_PAV;
% FALSE_POSITIVES=N_CALLS_PBSV_PAV-TRUE_POSITIVES;
% TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
% FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
% for i=[1:N_PANELS]
%     subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'x');
% endfor
%
% # SNIFFLES+PAV, all types.
% TRUE_POSITIVES=N_TP_SNIFFLES_PAV;
% FALSE_POSITIVES=N_CALLS_SNIFFLES_PAV-TRUE_POSITIVES;
% TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
% FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
% for i=[1:N_PANELS]
%     subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'s');
% endfor
%
% # PBSV, all types.
% TRUE_POSITIVES=N_TP_PBSV;
% FALSE_POSITIVES=N_CALLS_PBSV-TRUE_POSITIVES;
% TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
% FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
% for i=[1:N_PANELS]
%     subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'^');
% endfor
%
% # SNIFFLES, all types.
% TRUE_POSITIVES=N_TP_SNIFFLES;
% FALSE_POSITIVES=N_CALLS_SNIFFLES-TRUE_POSITIVES;
% TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
% FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
% for i=[1:N_PANELS]
%     subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'>');
% endfor
%
% # PAV, all types.
% TRUE_POSITIVES=N_TP_PAV;
% FALSE_POSITIVES=N_CALLS_PAV-TRUE_POSITIVES;
% TRUE_NEGATIVES=N_CALLS_ALL-N_TP-FALSE_POSITIVES;
% FALSE_NEGATIVES=N_TP-TRUE_POSITIVES;
% for i=[1:N_PANELS]
%     subplot(2,4,i); plot(FALSE_POSITIVES(i)/(FALSE_POSITIVES(i)+TRUE_NEGATIVES(i)),TRUE_POSITIVES(i)/(TRUE_POSITIVES(i)+FALSE_NEGATIVES(i)),'v');
% endfor