import java.io.*;


/**
 * Given a VCF of unsupported calls from $truvari collapse$, the program reports
 * the number of calls supported by each caller, using the SUPP annotation in 
 * the GT columns.
 */
public class AnalyzeFPs {
    
    /**
     * Remark: the program prints all AF values to STDOUT.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final int LENGTH_MIN = Integer.parseInt(args[1]);
        final int LENGTH_MAX = Integer.parseInt(args[2]);
        
        final char COMMENT = '#';
        final String GT_SEPARATOR = ":";
        final String SUPP_STR = "SUPP";
        final String GT_STR = "GT";
        final int MASK_PBSV = 1;
        final int MASK_SNIFFLES = 2;
        final int MASK_PAV = 4;
        
        int i, p;
        int length, supp, suppIndex, gtIndex;
        int nCalls, nCalls_noCaller, nCalls_sniffles, nCalls_pbsv, nCalls_pav, nCalls_sniffles_pbsv, nCalls_sniffles_pav, nCalls_pbsv_pav, nCalls_sniffles_pbsv_pav, nCalls_geq_two_callers;
        int nTotalCalls_sniffles, nTotalCalls_pbsv, nTotalCalls_pav;
        double af, nSamples;
        String str, field;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        nCalls=0; nCalls_noCaller=0; nCalls_sniffles=0; nCalls_pbsv=0; nCalls_pav=0; nCalls_sniffles_pbsv=0; nCalls_sniffles_pav=0; nCalls_pbsv_pav=0; nCalls_sniffles_pbsv_pav=0; nCalls_geq_two_callers=0;
        nTotalCalls_sniffles=0; nTotalCalls_pbsv=0; nTotalCalls_pav=0;
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                str=br.readLine();
                continue;
            }
            if (nCalls%100000==0) System.err.println("Processed "+nCalls+" calls");
            tokens=str.split("\t");
            
            field=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (field!=null) length=Integer.parseInt(field);
            else length=tokens[3].length()-tokens[4].length();
            if (length<0) length=-length;
            if (length<LENGTH_MIN || length>LENGTH_MAX) { str=br.readLine(); continue; }
            
            nCalls++;
            supp=0;
            
            // Finding the SUPP field
            suppIndex=-1; gtIndex=-1;
            tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(SUPP_STR)) suppIndex=i;
                else if (tokensPrime[i].equalsIgnoreCase(GT_STR)) gtIndex=i;
            }
            if (suppIndex==-1) {
                System.err.println("ERROR: No SUPP field in record: "+str);
                System.exit(1);
            }
            if (gtIndex==-1) {
                System.err.println("ERROR: No GT field in record: "+str);
                System.exit(1);
            }
            
            // Computing AF, and the union of all SUPP fields.
            supp=0; af=0; nSamples=0;
            for (i=9; i<tokens.length; i++) {
                tokensPrime=tokens[i].split(GT_SEPARATOR);
                if (tokensPrime[suppIndex].charAt(0)!='.') supp|=Integer.parseInt(tokensPrime[suppIndex]);
                af+=(tokensPrime[gtIndex].charAt(0)=='1'?1:0)+(tokensPrime[gtIndex].charAt(2)=='1'?1:0);
                nSamples+=((tokensPrime[gtIndex].charAt(0)=='1')||(tokensPrime[gtIndex].charAt(2)=='1'))?1:0;
            }
            if (supp==MASK_SNIFFLES) nCalls_sniffles++;
            else if (supp==MASK_PBSV) nCalls_pbsv++;
            else if (supp==MASK_PAV) nCalls_pav++;
            else if (supp==(MASK_PBSV|MASK_SNIFFLES|MASK_PAV)) nCalls_sniffles_pbsv_pav++;
            else if ((supp&MASK_SNIFFLES)!=0 && (supp&MASK_PBSV)!=0) nCalls_sniffles_pbsv++;
            else if ((supp&MASK_SNIFFLES)!=0 && (supp&MASK_PAV)!=0) nCalls_sniffles_pav++;
            else if ((supp&MASK_PBSV)!=0 && (supp&MASK_PAV)!=0) nCalls_pbsv_pav++;
            else nCalls_noCaller++;
            if ((supp&MASK_SNIFFLES)!=0) nTotalCalls_sniffles++;
            if ((supp&MASK_PBSV)!=0) nTotalCalls_pbsv++;
            if ((supp&MASK_PAV)!=0) nTotalCalls_pav++;
            if (Integer.bitCount(supp)>=2) nCalls_geq_two_callers++;
            System.out.println(""+af+","+(af/((tokens.length-9)*2))+","+nSamples);
            str=br.readLine();
        }
        br.close();

        System.err.println("sniffles only: "+((100.0*nCalls_sniffles)/nCalls)+"%");
        System.err.println("pbsv only: "+((100.0*nCalls_pbsv)/nCalls)+"%");
        System.err.println("pav only: "+((100.0*nCalls_pav)/nCalls)+"%");
        System.err.println("sniffles and pbsv: "+((100.0*nCalls_sniffles_pbsv)/nCalls)+"%");
        System.err.println("sniffles and pav: "+((100.0*nCalls_sniffles_pav)/nCalls)+"%");
        System.err.println("all callers: "+((100.0*nCalls_sniffles_pbsv_pav)/nCalls)+"%");
        System.err.println("no caller: "+((100.0*nCalls_noCaller)/nCalls)+"%");
        System.err.println(">=2 callers: "+((100.0*nCalls_geq_two_callers)/nCalls)+"%");
        System.err.println();
        System.err.println("All calls supported by sniffles: "+nTotalCalls_sniffles);
        System.err.println("All calls supported by pbsv: "+nTotalCalls_pbsv);
        System.err.println("All calls supported by pav: "+nTotalCalls_pav);
        System.err.println("All calls supported by >=2 callers: "+nCalls_geq_two_callers);
        System.err.println("All calls: "+nCalls);
    }
    
}