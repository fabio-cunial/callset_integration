import java.io.*;


/**
 * Testing GQ/DV correlation
 */
public class GQcurve8 {
    
    private static final char COMMENT = '#';
    private static final String GT_SEPARATOR = ":";
    private static final String GQ_STR = "GQ";
    private static final String DV_STR = "DV";
    
    private static int N_TRUE_CALLS, N_TRUE_DEL, N_TRUE_INS;
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String CALLER_VCF = args[0];  // All calls from the caller
        final String CALLER_TP_VCF = args[1];  // True positives from the caller
        final String TRUTH_VCF = args[2];
        final boolean GEQ = Integer.parseInt(args[3])==1;  // Unused
        final String OUTPUT_DIR = args[4];
        
        final int LENGTH_MIN = 50;
        final int LENGTH_MAX = 50000;
        
        int i;
        BufferedWriter bwAll, bwDel, bwIns, bwAll_tp, bwDel_tp, bwIns_tp;
        
        bwAll = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/gq_dv.txt"));
        bwDel = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/gq_dv_del.txt"));
        bwIns = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/gq_dv_ins.txt"));
        countCallsAt(CALLER_VCF,LENGTH_MIN,LENGTH_MAX,bwAll,bwDel,bwIns);
        bwAll.close(); bwDel.close(); bwIns.close();
        
        bwAll_tp = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/gq_dv_tp.txt"));
        bwDel_tp = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/gq_dv_del_tp.txt"));
        bwIns_tp = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/gq_dv_ins_tp.txt"));
        countCallsAt(CALLER_TP_VCF,LENGTH_MIN,LENGTH_MAX,bwAll_tp,bwDel_tp,bwIns_tp);
        bwAll_tp.close(); bwDel_tp.close(); bwIns_tp.close();
    }
    
    
    /**
     * For every $x$, counts the number of calls with $DR+DV>=x$ (if $gqGeq=
     * true$) or with $DR+DV<=x$ (if $gqGeq=false$).
     *
     * @param *_histogram stores the number of calls with DR+DV exactly equal to
     * a specific value (set to NULL to discard).
     */
    private static final void countCallsAt(String vcfFile, int lengthMin, int lengthMax, BufferedWriter bwAll, BufferedWriter bwDel, BufferedWriter bwIns) throws IOException {
        int i;
        int out, gqIndex, dvIndex, type, length;
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
            
            gqIndex=-1; dvIndex=-1;
            tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(GQ_STR)) gqIndex=i;
                else if (tokensPrime[i].equalsIgnoreCase(DV_STR)) dvIndex=i;
            }
            if (gqIndex==-1 || dvIndex==-1) { 
                System.err.println("This call has no GQ or DV: "+str);
                str=br.readLine();
                continue;
            }
            tokensPrime=tokens[9].split(GT_SEPARATOR);
            
            // All calls
            bwAll.write(tokensPrime[gqIndex]+","+tokensPrime[dvIndex]+"\n");
            
            // By SV type
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==VCFconstants.TYPE_DELETION) bwDel.write(tokensPrime[gqIndex]+","+tokensPrime[dvIndex]+"\n");
            else if (type==VCFconstants.TYPE_INSERTION) bwIns.write(tokensPrime[gqIndex]+","+tokensPrime[dvIndex]+"\n");
            str=br.readLine();
        }
    }
    
}