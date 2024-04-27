import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class ROCcurve {
    
    private static final char COMMENT = '#';
    private static char MISSING_CHR = '.';
    private static final String GT_SEPARATOR = ":";
    private static final String SCORE_STR = "HAPESTRY_READS_MAX";
    private static final String AD_STR = "AD";
    private static final String AD_SEPARATOR = ",";
    private static String TAG_STR;
    
    private static int[] N_TRUE_CALLS, N_TRUE_DEL, N_TRUE_INS;
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        TAG_STR=args[0];
        final boolean TAG_GEQ = Integer.parseInt(args[1])==1;
        final int LENGTH_MIN = Integer.parseInt(args[2]);
        final int LENGTH_MAX = Integer.parseInt(args[3]);
        final String CALLER_VCF = args[4];  // All calls from the caller
        final String CALLER_TP_VCF = args[5];  // True positives from the caller
        final String TRUTH_VCF = args[6];
        final String OUTPUT_DIR = args[7];
        final String TR_STATUS = args[8];
        final String GENOTYPER = args[9];
        final String SVLEN_BINS = args[10];  // Comma-separated, sorted.
        
        int i, j;
        int maxTAG, truePositives, falsePositives, trueNegatives, falseNegatives;
        int[] svlengths;
        String[] tokens;
        int[][] nCalls, nDel, nIns, nTP, nTP_del, nTP_ins, nCalls_histogram, nDel_histogram, nIns_histogram;
        BufferedWriter bw;
        
        tokens=SVLEN_BINS.split(",");
        svlengths = new int[tokens.length];
        for (i=0; i<tokens.length; i++) svlengths[i]=Integer.parseInt(tokens[i]);
        
        N_TRUE_CALLS = new int[2+svlengths.length];
        N_TRUE_DEL = new int[2+svlengths.length];
        N_TRUE_INS = new int[2+svlengths.length];
        getNTrueCalls(TRUTH_VCF,LENGTH_MIN,LENGTH_MAX,svlengths);
        maxTAG=getMaxTAG(CALLER_VCF);
        System.err.println("maxTAG="+maxTAG);
        nCalls = new int[2+svlengths.length][maxTAG+1];
        nDel = new int[2+svlengths.length][maxTAG+1];
        nIns = new int[2+svlengths.length][maxTAG+1];
        nCalls_histogram = new int[2+svlengths.length][maxTAG+1];
        nDel_histogram = new int[2+svlengths.length][maxTAG+1];
        nIns_histogram = new int[2+svlengths.length][maxTAG+1];
        countCallsAtTAG(CALLER_VCF,TAG_GEQ,LENGTH_MIN,LENGTH_MAX,svlengths,nCalls,nDel,nIns,nCalls_histogram,nDel_histogram,nIns_histogram);
        nTP = new int[2+svlengths.length][maxTAG+1];
        nTP_del = new int[2+svlengths.length][maxTAG+1];
        nTP_ins = new int[2+svlengths.length][maxTAG+1];
        countCallsAtTAG(CALLER_TP_VCF,TAG_GEQ,LENGTH_MIN,LENGTH_MAX,svlengths,nTP,nTP_del,nTP_ins,null,null,null);
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_roc_calls_"+TAG_STR+"_"+(TAG_GEQ?"geq":"leq")+".log"));
        for (i=0; i<=maxTAG; i++) {
            truePositives=nTP[0][i];
            falsePositives=nCalls[0][i]-truePositives;
            trueNegatives=nCalls[0][TAG_GEQ?0:nCalls.length-1]-nTP[0][TAG_GEQ?0:nCalls.length-1]-falsePositives;
            falseNegatives=nTP[0][TAG_GEQ?0:nCalls.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        for (j=0; j<svlengths.length; j++) {
            bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_"+svlengths[j]+"_roc_calls_"+TAG_STR+"_"+(TAG_GEQ?"geq":"leq")+".log"));
            for (i=0; i<=maxTAG; i++) {
                truePositives=nTP[j+1][i];
                falsePositives=nCalls[j+1][i]-truePositives;
                trueNegatives=nCalls[j+1][TAG_GEQ?0:nCalls.length-1]-nTP[j+1][TAG_GEQ?0:nCalls.length-1]-falsePositives;
                falseNegatives=nTP[j+1][TAG_GEQ?0:nCalls.length-1]-truePositives;
                bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
            }
            bw.close();
        }
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_roc_del_"+TAG_STR+"_"+(TAG_GEQ?"geq":"leq")+".log"));
        for (i=0; i<=maxTAG; i++) {
            truePositives=nTP_del[0][i];
            falsePositives=nDel[0][i]-truePositives;
            trueNegatives=nDel[0][TAG_GEQ?0:nDel.length-1]-nTP_del[0][TAG_GEQ?0:nTP_del.length-1]-falsePositives;
            falseNegatives=nTP_del[0][TAG_GEQ?0:nTP_del.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        for (j=0; j<svlengths.length; j++) {
            bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_"+svlengths[j]+"_roc_del_"+TAG_STR+"_"+(TAG_GEQ?"geq":"leq")+".log"));
            for (i=0; i<=maxTAG; i++) {
                truePositives=nTP_del[j+1][i];
                falsePositives=nDel[j+1][i]-truePositives;
                trueNegatives=nDel[j+1][TAG_GEQ?0:nDel.length-1]-nTP_del[j+1][TAG_GEQ?0:nTP_del.length-1]-falsePositives;
                falseNegatives=nTP_del[j+1][TAG_GEQ?0:nTP_del.length-1]-truePositives;
                bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
            }
            bw.close();
        }
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_roc_ins_"+TAG_STR+"_"+(TAG_GEQ?"geq":"leq")+".log"));
        for (i=0; i<=maxTAG; i++) {
            truePositives=nTP_ins[0][i];
            falsePositives=nIns[0][i]-truePositives;
            trueNegatives=nIns[0][TAG_GEQ?0:nIns.length-1]-nTP_ins[0][TAG_GEQ?0:nTP_ins.length-1]-falsePositives;
            falseNegatives=nTP_ins[0][TAG_GEQ?0:nTP_ins.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        for (j=0; j<svlengths.length; j++) {
            bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_"+svlengths[j]+"_roc_ins_"+TAG_STR+"_"+(TAG_GEQ?"geq":"leq")+".log"));
            for (i=0; i<=maxTAG; i++) {
                truePositives=nTP_ins[j+1][i];
                falsePositives=nIns[j+1][i]-truePositives;
                trueNegatives=nIns[j+1][TAG_GEQ?0:nIns.length-1]-nTP_ins[j+1][TAG_GEQ?0:nTP_ins.length-1]-falsePositives;
                falseNegatives=nTP_ins[j+1][TAG_GEQ?0:nTP_ins.length-1]-truePositives;
                bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
            }
            bw.close();
        }
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_histogram_calls_"+TAG_STR+".log"));
        for (i=0; i<=maxTAG; i++) bw.write(i+","+nCalls_histogram[0][i]+"\n");
        bw.close();
        for (j=0; j<svlengths.length; j++) {
            bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_"+svlengths[j]+"_histogram_calls_"+TAG_STR+".log"));
            for (i=0; i<=maxTAG; i++) bw.write(i+","+nCalls_histogram[j+1][i]+"\n");
            bw.close();
        }
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_histogram_del_"+TAG_STR+".log"));
        for (i=0; i<=maxTAG; i++) bw.write(i+","+nDel_histogram[0][i]+"\n");
        bw.close();
        for (j=0; j<svlengths.length; j++) {
            bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_"+svlengths[j]+"_histogram_del_"+TAG_STR+".log"));
            for (i=0; i<=maxTAG; i++) bw.write(i+","+nDel_histogram[j+1][i]+"\n");
            bw.close();
        }
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_histogram_ins_"+TAG_STR+".log"));
        for (i=0; i<=maxTAG; i++) bw.write(i+","+nIns_histogram[0][i]+"\n");
        bw.close();
        for (j=0; j<svlengths.length; j++) {
            bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_"+svlengths[j]+"_histogram_ins_"+TAG_STR+".log"));
            for (i=0; i<=maxTAG; i++) bw.write(i+","+nIns_histogram[j+1][i]+"\n");
            bw.close();
        }
    }
    
    
    /**
     * Sets global arrays $N_TRUE_CALLS,N_TRUE_DEL,N_TRUE_INS$.
     */
    private static final void getNTrueCalls(String vcfFile, int lengthMin, int lengthMax, int[] svlengths) throws IOException {
        int i;
        int type, length, binID;
        String str, field;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(vcfFile));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) { str=br.readLine(); continue; }
            tokens=str.split("\t");
            
            field=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (field!=null) length=Integer.parseInt(field);
            else length=tokens[3].length()-tokens[4].length();
            if (length<0) length=-length;
            if (length<lengthMin || length>lengthMax) { str=br.readLine(); continue; }
            binID=Arrays.binarySearch(svlengths,length);
            if (binID<0) binID=-1-binID;
            
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==-1) {
                if (tokens[3].length()>tokens[4].length()) type=VCFconstants.TYPE_DELETION;
                else if (tokens[3].length()<tokens[4].length()) type=VCFconstants.TYPE_INSERTION;
            }
            
            N_TRUE_CALLS[0]++; N_TRUE_CALLS[1+binID]++;
            if (type==VCFconstants.TYPE_DELETION) { N_TRUE_DEL[0]++; N_TRUE_DEL[1+binID]++; }
            else if (type==VCFconstants.TYPE_INSERTION) { N_TRUE_INS[0]++; N_TRUE_INS[1+binID]++; }
            
            str=br.readLine();
        }
        
        System.err.println("N_TRUE_CALLS:");
        System.err.println("-1,"+N_TRUE_CALLS[0]);
        for (i=0; i<svlengths.length; i++) System.err.println(svlengths[i]+","+N_TRUE_CALLS[i+1]);
        System.err.println("N_TRUE_DEL:");
        System.err.println("-1,"+N_TRUE_DEL[0]);
        for (i=0; i<svlengths.length; i++) System.err.println(svlengths[i]+","+N_TRUE_DEL[i+1]);
        System.err.println("N_TRUE_INS:");
        System.err.println("-1,"+N_TRUE_INS[0]);
        for (i=0; i<svlengths.length; i++) System.err.println(svlengths[i]+","+N_TRUE_INS[i+1]);
    }
    
    
    /**
     * @return the max TAG value over all records in $vcfFile$.
     */
    private static final int getMaxTAG(String vcfFile) throws IOException {
        int i, p;
        int out, tagIndex, value;
        String str, field;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        out=0;
        br = new BufferedReader(new FileReader(vcfFile));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) { str=br.readLine(); continue; }
            tokens=str.split("\t");
            if (TAG_STR.equals(SCORE_STR)) {
                field=VCFconstants.getField(tokens[7],SCORE_STR);
                if (field==null) { str=br.readLine(); continue; }
                value=(int)(1000*Double.parseDouble(field));
            }
            else {
                tagIndex=-1;
                tokensPrime=tokens[8].split(GT_SEPARATOR);
                for (i=0; i<tokensPrime.length; i++) {
                    if (tokensPrime[i].equalsIgnoreCase(TAG_STR)) { tagIndex=i; break; }
                }
                if (tagIndex==-1) { str=br.readLine(); continue; }
                tokensPrime=tokens[9].split(GT_SEPARATOR);
                if (tokensPrime[tagIndex].charAt(0)==MISSING_CHR) value=0;
                else if (TAG_STR.equals(AD_STR)) {
                    p=tokensPrime[tagIndex].indexOf(AD_SEPARATOR);
                    value=(int)Double.parseDouble(tokensPrime[tagIndex].substring(p+1));
                }
                else value=Integer.parseInt(tokensPrime[tagIndex]);
            }
            if (value>out) out=value;
            str=br.readLine();
        }
        return out;
    }
    
    
    /**
     * For every $x$, counts the number of calls with $TAG>=x$ (if 
     * $valueGeq=true$) or with $TAG<=x$ (if $valueGeq=false$).
     *
     * @param *_histogram stores the number of calls with TAG exactly equal to a
     * specific value (set to NULL to discard).
     */
    private static final void countCallsAtTAG(String vcfFile, boolean valueGeq, int lengthMin, int lengthMax, int[] svlengths, int[][] nCalls, int[][] nDel, int[][] nIns, int[][] nCalls_histogram, int[][] nDel_histogram, int[][] nIns_histogram) throws IOException {
        int i, p;
        int out, tagIndex, value, type, length, binID;
        String str, field;
        BufferedReader br;
        String[] tokens, tokensPrime;
                  
        br = new BufferedReader(new FileReader(vcfFile));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) { str=br.readLine(); continue; }
            tokens=str.split("\t");
            
            field=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (field!=null) length=Integer.parseInt(field);
            else length=tokens[3].length()-tokens[4].length();
            if (length<0) length=-length;
            if (length<lengthMin || length>lengthMax) { str=br.readLine(); continue; }
            binID=Arrays.binarySearch(svlengths,length);
            if (binID<0) binID=-1-binID;
            
            if (TAG_STR.equals(SCORE_STR)) {
                field=VCFconstants.getField(tokens[7],SCORE_STR);
                if (field==null) { 
                    System.err.println("This call has no TAG: "+str);
                    str=br.readLine();
                    continue;
                }
                value=(int)(1000*Double.parseDouble(field));
            }
            else {
                tagIndex=-1;
                tokensPrime=tokens[8].split(GT_SEPARATOR);
                for (i=0; i<tokensPrime.length; i++) {
                    if (tokensPrime[i].equalsIgnoreCase(TAG_STR)) { tagIndex=i; break; }
                }
                if (tagIndex==-1) { 
                    System.err.println("This call has no TAG: "+str);
                    str=br.readLine();
                    continue;               
                }
                tokensPrime=tokens[9].split(GT_SEPARATOR);
                if (tokensPrime[tagIndex].charAt(0)==MISSING_CHR) value=0;
                else if (TAG_STR.equals(AD_STR)) {
                    p=tokensPrime[tagIndex].indexOf(AD_SEPARATOR);
                    value=(int)Double.parseDouble(tokensPrime[tagIndex].substring(p+1));
                }
                else value=Integer.parseInt(tokensPrime[tagIndex]);
            }
            
            // All calls
            if (nCalls_histogram!=null) { nCalls_histogram[0][value]++; nCalls_histogram[1+binID][value]++; }
            if (valueGeq) {
                for (i=0; i<=value; i++) { nCalls[0][i]++; nCalls[1+binID][i]++; }
            }
            else {
                for (i=value; i<nCalls.length; i++) { nCalls[0][i]++; nCalls[1+binID][i]++; }
            }
            
            // By SV type
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==VCFconstants.TYPE_DELETION) {
                if (nDel_histogram!=null) { nDel_histogram[0][value]++; nDel_histogram[1+binID][value]++; }
                if (valueGeq) {
                    for (i=0; i<=value; i++) { nDel[0][i]++; nDel[1+binID][i]++; }
                }
                else {
                    for (i=value; i<nDel.length; i++) { nDel[0][i]++; nDel[1+binID][i]++; }
                }
            }
            else if (type==VCFconstants.TYPE_INSERTION) {
                if (nIns_histogram!=null) { nIns_histogram[0][value]++; nIns_histogram[1+binID][value]++; }
                if (valueGeq) {
                    for (i=0; i<=value; i++) { nIns[0][i]++; nIns[1+binID][i]++; }
                }
                else {
                    for (i=value; i<nIns.length; i++) { nIns[0][i]++; nIns[1+binID][i]++; }
                }
            }
            str=br.readLine();
        }
    }
    
}