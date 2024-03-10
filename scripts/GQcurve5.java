import java.io.*;


/**
 * 
 */
public class GQcurve5 {
    
    private static final char COMMENT = '#';
    private static final String GT_SEPARATOR = ":";
    private static final String AD_SEPARATOR = ",";
    private static String AD_STR;
    
    private static int N_TRUE_CALLS, N_TRUE_DEL, N_TRUE_INS;
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        AD_STR=args[0];
        final boolean GEQ = Integer.parseInt(args[1])==1;
        final int LENGTH_MIN = Integer.parseInt(args[2]);
        final int LENGTH_MAX = Integer.parseInt(args[3]);
        final String CALLER_VCF = args[4];  // All calls from the caller
        final String CALLER_TP_VCF = args[5];  // True positives from the caller
        final String TRUTH_VCF = args[6];        
        final String OUTPUT_DIR = args[7];
        final String TR_STATUS = args[8];
        final String GENOTYPER = args[9];
        
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
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_roc_calls_"+AD_STR+"_"+(GEQ?"geq":"leq")+".txt"));
        for (i=0; i<=max; i++) {
            truePositives=nTP[i];
            falsePositives=nCalls[i]-truePositives;
            trueNegatives=nCalls[GEQ?0:nCalls.length-1]-nTP[GEQ?0:nCalls.length-1]-falsePositives;
            falseNegatives=nTP[GEQ?0:nCalls.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_roc_del_"+AD_STR+"_"+(GEQ?"geq":"leq")+".txt"));
        for (i=0; i<=max; i++) {
            truePositives=nTP_del[i];
            falsePositives=nDel[i]-truePositives;
            trueNegatives=nDel[GEQ?0:nDel.length-1]-nTP_del[GEQ?0:nTP_del.length-1]-falsePositives;
            falseNegatives=nTP_del[GEQ?0:nTP_del.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_roc_ins_"+AD_STR+"_"+(GEQ?"geq":"leq")+".txt"));
        for (i=0; i<=max; i++) {
            truePositives=nTP_ins[i];
            falsePositives=nIns[i]-truePositives;
            trueNegatives=nIns[GEQ?0:nIns.length-1]-nTP_ins[GEQ?0:nTP_ins.length-1]-falsePositives;
            falseNegatives=nTP_ins[GEQ?0:nTP_ins.length-1]-truePositives;
            bw.write(i+","+(((double)falsePositives)/(falsePositives+trueNegatives))+","+(((double)truePositives)/(truePositives+falseNegatives))+"\n");
        }
        bw.close();
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_histogram_calls_"+AD_STR+".txt"));
        for (i=0; i<=max; i++) bw.write(i+","+nCalls_histogram[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_histogram_del_"+AD_STR+".txt"));
        for (i=0; i<=max; i++) bw.write(i+","+nDel_histogram[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+TR_STATUS+"_"+GENOTYPER+"_histogram_ins_"+AD_STR+".txt"));
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
     * @return the max value over all records in $vcfFile$.
     */
    private static final int getMax(String vcfFile) throws IOException {
        int i, p;
        int out, index, value;
        String str, strPrime;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        out=0;
        br = new BufferedReader(new FileReader(vcfFile));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) { str=br.readLine(); continue; }
            tokens=str.split("\t");
            index=-1;
            tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(AD_STR)) { index=i; break; }
            }
            if (index==-1) { str=br.readLine(); continue; }
            tokensPrime=tokens[9].split(GT_SEPARATOR);
            strPrime=tokensPrime[index];
            p=strPrime.indexOf(AD_SEPARATOR);
            value=strPrime.substring(p+1).charAt(0)=='.'?0:(int)Double.parseDouble(strPrime.substring(p+1));
            if (value>out) out=value;
            str=br.readLine();
        }
        return out;
    }
    
    
    /**
     * For every $x$, counts the number of calls with value>=x (if $geq=true$)
     * or with value<=x (if $geq=false$).
     *
     * @param *_histogram stores the number of calls with value exactly equal to
     * a specific value (set to NULL to discard).
     */
    private static final void countCallsAt(String vcfFile, boolean geq, int lengthMin, int lengthMax, int[] nCalls, int[] nDel, int[] nIns, int[] nCalls_histogram, int[] nDel_histogram, int[] nIns_histogram) throws IOException {
        int i, p;
        int out, index, value, type, length;
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
            
            index=-1;
            tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(AD_STR)) { index=i; break; }
            }
            if (index==-1) { 
                System.err.println("This call has no AD: "+str);
                str=br.readLine();
                continue;
            }
            tokensPrime=tokens[9].split(GT_SEPARATOR);
            field=tokensPrime[index];
            p=field.indexOf(AD_SEPARATOR);
            value=field.substring(p+1).charAt(0)=='.'?0:(int)Double.parseDouble(field.substring(p+1));
            
            // All calls
            if (nCalls_histogram!=null) nCalls_histogram[value]++;
            if (geq) {
                for (i=0; i<=value; i++) nCalls[i]++;
            }
            else {
                for (i=value; i<nCalls.length; i++) nCalls[i]++;
            }
            
            // By SV type
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==VCFconstants.TYPE_DELETION) {
                if (nDel_histogram!=null) nDel_histogram[value]++;
                if (geq) {
                    for (i=0; i<=value; i++) nDel[i]++;
                }
                else {
                    for (i=value; i<nDel.length; i++) nDel[i]++;
                }
            }
            else if (type==VCFconstants.TYPE_INSERTION) {
                if (nIns_histogram!=null) nIns_histogram[value]++;
                if (geq) {
                    for (i=0; i<=value; i++) nIns[i]++;
                }
                else {
                    for (i=value; i<nIns.length; i++) nIns[i]++;
                }
            }
            str=br.readLine();
        }
    }
    
}