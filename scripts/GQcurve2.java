import java.io.*;


/**
 * 
 */
public class GQcurve2 {
    
    private static final char COMMENT = '#';
    private static final String GT_SEPARATOR = ":";
    private static final String DR_STR = "DR";
    private static final String DV_STR = "DV";
    
    private static int N_TRUE_CALLS, N_TRUE_DEL, N_TRUE_INS;
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final boolean GEQ = Integer.parseInt(args[0])==1;
        final int LENGTH_MIN = Integer.parseInt(args[1]);
        final int LENGTH_MAX = Integer.parseInt(args[2]);
        final String CALLER_VCF = args[3];  // All calls from the caller
        final String CALLER_TP_VCF = args[4];  // True positives from the caller
        final String TRUTH_VCF = args[5];
        final String OUTPUT_DIR = args[6];
        
        int i;
        int max, truePositives, falsePositives, trueNegatives, falseNegatives;
        int[] nCalls, nDel, nIns, nTP, nTP_del, nTP_ins, nCalls_histogram, nDel_histogram, nIns_histogram;
        BufferedWriter bw;
        
        getNTrueCalls(TRUTH_VCF,LENGTH_MIN,LENGTH_MAX);
        max=getMax(CALLER_VCF);
        System.err.println("max="+max);
        nCalls = new int[max+1];
        nDel = new int[max+1];
        nIns = new int[max+1];
        nCalls_histogram = new int[max+1];
        nDel_histogram = new int[max+1];
        nIns_histogram = new int[max+1];
        countCallsAt(CALLER_VCF,GEQ,LENGTH_MIN,LENGTH_MAX,nCalls,nDel,nIns,nCalls_histogram,nDel_histogram,nIns_histogram);
        nTP = new int[max+1];
        nTP_del = new int[max+1];
        nTP_ins = new int[max+1];
        countCallsAt(CALLER_TP_VCF,GEQ,LENGTH_MIN,LENGTH_MAX,nTP,nTP_del,nTP_ins,null,null,null);
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/roc_calls_DRDV_"+(GEQ?"geq":"leq")+".txt"));
        for (i=0; i<=max; i++) {
            truePositives=nTP[i];
            falsePositives=nCalls[i]-truePositives;
            trueNegatives=nCalls[GEQ?0:nCalls.length-1]-nTP[GEQ?0:nCalls.length-1]-falsePositives;
            falseNegatives=nTP[GEQ?0:nCalls.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/roc_del_DRDV_"+(GEQ?"geq":"leq")+".txt"));
        for (i=0; i<=max; i++) {
            truePositives=nTP_del[i];
            falsePositives=nDel[i]-truePositives;
            trueNegatives=nDel[GEQ?0:nDel.length-1]-nTP_del[GEQ?0:nTP_del.length-1]-falsePositives;
            falseNegatives=nTP_del[GEQ?0:nTP_del.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/roc_ins_DRDV_"+(GEQ?"geq":"leq")+".txt"));
        for (i=0; i<=max; i++) {
            truePositives=nTP_ins[i];
            falsePositives=nIns[i]-truePositives;
            trueNegatives=nIns[GEQ?0:nIns.length-1]-nTP_ins[GEQ?0:nTP_ins.length-1]-falsePositives;
            falseNegatives=nTP_ins[GEQ?0:nTP_ins.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/histogram_calls_DRDV.txt"));
        for (i=0; i<=max; i++) bw.write(i+","+nCalls_histogram[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/histogram_del_DRDV.txt"));
        for (i=0; i<=max; i++) bw.write(i+","+nDel_histogram[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/histogram_ins_DRDV.txt"));
        for (i=0; i<=max; i++) bw.write(i+","+nIns_histogram[i]+"\n");
        bw.close();
    }
    
    
    /**
     * Sets global variables $N_TRUE_CALLS,N_TRUE_DEL,N_TRUE_INS$.
     */
    private static final void getNTrueCalls(String vcfFile, int lengthMin, int lengthMax) throws IOException {
        int i;
        int type, length;
        String str, field;
        BufferedReader br;
        String[] tokens;
        
        N_TRUE_CALLS=0; N_TRUE_DEL=0; N_TRUE_INS=0;
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
            
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==-1) {
                if (tokens[3].length()>tokens[4].length()) type=VCFconstants.TYPE_DELETION;
                else if (tokens[3].length()<tokens[4].length()) type=VCFconstants.TYPE_INSERTION;
            }
            
            N_TRUE_CALLS++;
            if (type==VCFconstants.TYPE_DELETION) N_TRUE_DEL++;
            else if (type==VCFconstants.TYPE_INSERTION) N_TRUE_INS++;
            
            str=br.readLine();
        }
        System.err.println("N_TRUE_CALLS="+N_TRUE_CALLS+" N_TRUE_DEL="+N_TRUE_DEL+" N_TRUE_INS="+N_TRUE_INS);
    }
    
    
    /**
     * @return the max DR+DV value over all records in $vcfFile$.
     */
    private static final int getMax(String vcfFile) throws IOException {
        int i;
        int out, drIndex, dvIndex, sum;
        String str;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        out=0;
        br = new BufferedReader(new FileReader(vcfFile));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) { str=br.readLine(); continue; }
            tokens=str.split("\t");
            drIndex=-1; dvIndex=-1;
            tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(DR_STR)) drIndex=i;
                else if (tokensPrime[i].equalsIgnoreCase(DV_STR)) dvIndex=i;
            }
            if (drIndex==-1 || dvIndex==-1) { str=br.readLine(); continue; }
            tokensPrime=tokens[9].split(GT_SEPARATOR);
            sum=Integer.parseInt(tokensPrime[drIndex])+Integer.parseInt(tokensPrime[dvIndex]);
            if (sum>out) out=sum;
            str=br.readLine();
        }
        return out;
    }
    
    
    /**
     * For every $x$, counts the number of calls with $DR+DV>=x$ (if $gqGeq=
     * true$) or with $DR+DV<=x$ (if $gqGeq=false$).
     *
     * @param *_histogram stores the number of calls with DR+DV exactly equal to
     * a specific value (set to NULL to discard).
     */
    private static final void countCallsAt(String vcfFile, boolean gqGeq, int lengthMin, int lengthMax, int[] nCalls, int[] nDel, int[] nIns, int[] nCalls_histogram, int[] nDel_histogram, int[] nIns_histogram) throws IOException {
        int i;
        int out, drIndex, dvIndex, sum, type, length;
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
            
            drIndex=-1; dvIndex=-1;
            tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(DR_STR)) drIndex=i;
                else if (tokensPrime[i].equalsIgnoreCase(DV_STR)) dvIndex=i;
            }
            if (drIndex==-1 || dvIndex==-1) { 
                System.err.println("This call has no DR or DV: "+str);
                str=br.readLine();
                continue;
            }
            tokensPrime=tokens[9].split(GT_SEPARATOR);
            sum=Integer.parseInt(tokensPrime[drIndex])+Integer.parseInt(tokensPrime[dvIndex]);
            
            // All calls
            if (nCalls_histogram!=null) nCalls_histogram[sum]++;
            if (gqGeq) {
                for (i=0; i<=sum; i++) nCalls[i]++;
            }
            else {
                for (i=sum; i<nCalls.length; i++) nCalls[i]++;
            }
            
            // By SV type
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==VCFconstants.TYPE_DELETION) {
                if (nDel_histogram!=null) nDel_histogram[sum]++;
                if (gqGeq) {
                    for (i=0; i<=sum; i++) nDel[i]++;
                }
                else {
                    for (i=sum; i<nDel.length; i++) nDel[i]++;
                }
            }
            else if (type==VCFconstants.TYPE_INSERTION) {
                if (nIns_histogram!=null) nIns_histogram[sum]++;
                if (gqGeq) {
                    for (i=0; i<=sum; i++) nIns[i]++;
                }
                else {
                    for (i=sum; i<nIns.length; i++) nIns[i]++;
                }
            }
            str=br.readLine();
        }
    }
    
}