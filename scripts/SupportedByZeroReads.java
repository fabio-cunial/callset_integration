import java.util.Arrays;
import java.io.*;


/**
 * Given an intra-sample-merged VCF from $truvari collapse$, that was then
 * re-genotyped, the program counts calls supported by zero reads (e.g with 
 * DV=0). In particular, the program prints to STDOUT one line for each SVLEN
 * bin, with the following comma-separated values:
 *
 * sv_length_bin
 *
 * 0: nCalls
 * 1: nCalls supported exclusively by pbsv
 * 2: nCalls supported exclusively by sniffles
 * 3: nCalls supported exclusively by pav
 *
 * 4: nDel
 * 5: nDel supported exclusively by pbsv
 * 6: nDel supported exclusively by sniffles
 * 7: nDel supported exclusively by pav
 *
 * 8: nIns
 * 9: nIns supported exclusively by pbsv
 * 10: nIns supported exclusively by sniffles
 * 11: nIns supported exclusively by pav
 *
 * 12: nCalls supported by >1 caller
 * 13: nDel supported by >1 caller
 * 14: nIns supported by >1 caller
 *
 * 15: nCalls supported by zero callers
 * 16: nDel supported by zero callers
 * 17: nIns supported by zero callers
 *
 * The first line of the file contains totals over all length bins.
 */
public class SupportedByZeroReads {
    
    private static final char COMMENT = '#';
    private static final char MISSING_CHR = '.';
    private static final String GT_SEPARATOR = ":";
    private static final String AD_STR = "AD";
    private static final String AD_SEPARATOR = ",";
    private static final String SUPP_STR = "SUPP";
    private static final int MASK_PBSV = 1;
    private static final int MASK_SNIFFLES = 2;
    private static final int MASK_PAV = 4;
    
    
    public static void main(String[] args) throws IOException {
        final String TAG_STR = args[0];
        final String MERGED_VCF = args[1];
        final String REGENOTYPED_VCF = args[2];
        final int LENGTH_MIN = Integer.parseInt(args[3]);
        final int LENGTH_MAX = Integer.parseInt(args[4]);
        final String SVLEN_BINS = args[5];  // Comma-separated, sorted.
        
        int i, j, p;
        int length, supp, suppIndex, value, tagIndex, type, binID;
        int nRows;
        String str1, str2, field;
        BufferedReader br1, br2;
        String[] tokens, tokensPrime;
        int[] svlengths;
        int[][] out;
        
        tokens=SVLEN_BINS.split(",");
        svlengths = new int[tokens.length];
        for (i=0; i<tokens.length; i++) svlengths[i]=Integer.parseInt(tokens[i]);
        out = new int[2+svlengths.length][18];
        br1 = new BufferedReader(new FileReader(REGENOTYPED_VCF));
        br2 = new BufferedReader(new FileReader(MERGED_VCF));
        str1=br1.readLine();  // Ignoring headers
        while (str1!=null && str1.charAt(0)==COMMENT) str1=br1.readLine();
        str2=br2.readLine();
        while (str2!=null && str2.charAt(0)==COMMENT) str2=br2.readLine();
        nRows=0;
        while (str1!=null) {
            if (nRows%100000==0) System.err.println("Processed "+nRows+" calls");
            nRows++;
            tokens=str1.split("\t");
            
            // Filtering by SVLEN
            field=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (field!=null) length=Integer.parseInt(field);
            else length=tokens[3].length()-tokens[4].length();
            if (length<0) length=-length;
            if (length<LENGTH_MIN || length>LENGTH_MAX) { str1=br1.readLine(); str2=br2.readLine(); continue; }
            binID=Arrays.binarySearch(svlengths,length);
            if (binID<0) binID=-1-binID;
            
            // Keeping only calls with zero supporting reads
            tagIndex=-1;
            tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(TAG_STR)) { tagIndex=i; break; }
            }
            if (tagIndex==-1) {
                System.err.println("This call has no "+TAG_STR+": "+str1);
                str1=br1.readLine(); str2=br2.readLine();
                continue;
            }
            tokensPrime=tokens[9].split(GT_SEPARATOR);
            if (tokensPrime[tagIndex].charAt(0)==MISSING_CHR) value=0;
            else if (TAG_STR.equals(AD_STR)) {
                p=tokensPrime[tagIndex].indexOf(AD_SEPARATOR);
                value=(int)Double.parseDouble(tokensPrime[tagIndex].substring(p+1));
            }
            else value=Integer.parseInt(tokensPrime[tagIndex]);
            if (value!=0) { str1=br1.readLine(); str2=br2.readLine(); continue; }
            out[0][0]++; out[1+binID][0]++;
            
            // Reading SVTYPE
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==VCFconstants.TYPE_DELETION) { out[0][4]++; out[1+binID][4]++; }
            else if (type==VCFconstants.TYPE_INSERTION) { out[0][8]++; out[1+binID][8]++; }
            
            // Reading SUPP
            tokens=str2.split("\t");
            suppIndex=-1; tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(SUPP_STR)) suppIndex=i;
            }
            if (suppIndex==-1) {
                System.err.println("ERROR: No SUPP field in record: "+str2);
                System.exit(1);
            }
            tokensPrime=tokens[9].split(GT_SEPARATOR);
            supp=tokensPrime[suppIndex].charAt(0)=='.'?0:Integer.parseInt(tokensPrime[suppIndex]);
            if (supp==0) {
                out[0][15]++; out[1+binID][15]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][16]++; out[1+binID][16]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][17]++; out[1+binID][17]++; }
            }
            else if (Integer.bitCount(supp)>=2) {
                out[0][12]++; out[1+binID][12]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][13]++; out[1+binID][13]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][14]++; out[1+binID][14]++; }
            }
            else if ((supp&MASK_PBSV)!=0 && (supp&MASK_SNIFFLES)==0 && (supp&MASK_PAV)==0) {
                out[0][1]++; out[1+binID][1]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][5]++; out[1+binID][5]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][9]++; out[1+binID][9]++; }
            }
            else if ((supp&MASK_PBSV)==0 && (supp&MASK_SNIFFLES)!=0 && (supp&MASK_PAV)==0) {
                out[0][2]++; out[1+binID][2]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][6]++; out[1+binID][6]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][10]++; out[1+binID][10]++; }
            }
            else if ((supp&MASK_PBSV)==0 && (supp&MASK_SNIFFLES)==0 && (supp&MASK_PAV)!=0) {
                out[0][3]++; out[1+binID][3]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][7]++; out[1+binID][7]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][11]++; out[1+binID][11]++; }
            }
            else System.err.println("ERROR: SUPP="+supp+" is not captured by any case.");
            
            // Next iteration
            str1=br1.readLine(); str2=br2.readLine();
        }
        br1.close(); br2.close();
        
        // Outputting
        System.out.print("-1,"+out[0][0]+"");
        for (j=1; j<out[0].length; j++) System.out.print(","+out[0][j]);
        System.out.println();
        for (i=0; i<svlengths.length; i++) {
            System.out.print(svlengths[i]+","+out[1+i][0]+"");
            for (j=1; j<out[1+i].length; j++) System.out.print(","+out[1+i][j]);
            System.out.println();
        }
    }
    
}