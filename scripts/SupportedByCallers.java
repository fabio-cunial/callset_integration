import java.util.Arrays;
import java.io.*;


/**
 * Given an intra-sample-merged VCF from $truvari collapse$, the program prints
 * to STDOUT one line for each SVLEN bin, with the following comma-separated 
 * values:
 *
 * sv_length_bin
 *
 * 0: nCalls
 * 1: nDel
 * 2: nIns
 *
 * 3: nCalls supported by >=2 callers
 * 4: nDel supported by >=2 callers
 * 5: nIns supported by >=2 callers
 *
 * 6: nCalls calls supported by 3 callers
 * 7: nDel calls supported by 3 callers
 * 8: nIns calls supported by 3 callers
 *
 * 9: nCalls supported by pbsv and sniffles
 * 10: nDel supported by pbsv and sniffles
 * 11: nIns supported by pbsv and sniffles
 *
 * 12: nCalls supported by pbsv and pav
 * 13: nDel supported by pbsv and pav
 * 14: nIns supported by pbsv and pav
 *
 * 15: nCalls supported by sniffles and pav
 * 16: nDel supported by sniffles and pav
 * 17: nIns supported by sniffles and pav
 *
 * 18: nCalls supported by pbsv (and possibly other callers)
 * 19: nDel supported by pbsv (and possibly other callers)
 * 20: nIns supported by pbsv (and possibly other callers)
 * 
 * 21: nCalls supported by sniffles (and possibly other callers)
 * 22: nDel supported by sniffles (and possibly other callers)
 * 23: nIns supported by sniffles (and possibly other callers)
 *
 * 24: nCalls supported by pav (and possibly other callers)
 * 25: nDel supported by pav (and possibly other callers)
 * 26: nIns supported by pav (and possibly other callers)
 */
public class SupportedByCallers {
    
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final int LENGTH_MIN = Integer.parseInt(args[1]);
        final int LENGTH_MAX = Integer.parseInt(args[2]);
        final String SVLEN_BINS = args[3];  // Comma-separated, sorted.
        
        final char COMMENT = '#';
        final String GT_SEPARATOR = ":";
        final String SUPP_STR = "SUPP";
        final String GT_STR = "GT";
        final int MASK_PBSV = 1;
        final int MASK_SNIFFLES = 2;
        final int MASK_PAV = 4;
        
        int i, j;
        int length, supp, suppIndex, type, binID, nRows;
        int[] svlengths;
        int[][] out;
        String str, field;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        tokens=SVLEN_BINS.split(",");
        svlengths = new int[tokens.length];
        for (i=0; i<tokens.length; i++) svlengths[i]=Integer.parseInt(tokens[i]);
        out = new int[2+svlengths.length][27];
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine(); nRows=0;
        while (str!=null) {
            if (str.charAt(0)==COMMENT) { str=br.readLine(); continue; }
            if (nRows%100000==0) System.err.println("Processed "+nRows+" calls");
            nRows++;
            tokens=str.split("\t");
            
            // Filtering by SVLEN
            field=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (field!=null) length=Integer.parseInt(field);
            else length=tokens[3].length()-tokens[4].length();
            if (length<0) length=-length;
            if (length<LENGTH_MIN || length>LENGTH_MAX) { str=br.readLine(); continue; }
            binID=Arrays.binarySearch(svlengths,length);
            if (binID<0) binID=-1-binID;
            
            out[0][0]++; out[1+binID][0]++;
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==VCFconstants.TYPE_DELETION) { out[0][1]++; out[1+binID][1]++; }
            else if (type==VCFconstants.TYPE_INSERTION) { out[0][2]++; out[1+binID][2]++; }
            
            suppIndex=-1; tokensPrime=tokens[8].split(GT_SEPARATOR);
            for (i=0; i<tokensPrime.length; i++) {
                if (tokensPrime[i].equalsIgnoreCase(SUPP_STR)) suppIndex=i;
            }
            if (suppIndex==-1) {
                System.err.println("ERROR: No SUPP field in record: "+str);
                System.exit(1);
            }
            tokensPrime=tokens[9].split(GT_SEPARATOR);
            supp=tokensPrime[suppIndex].charAt(0)=='.'?0:Integer.parseInt(tokensPrime[suppIndex]);
            if (Integer.bitCount(supp)==3) {
                out[0][6]++; out[1+binID][6]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][7]++; out[1+binID][7]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][8]++; out[1+binID][8]++; }
            }
            if (Integer.bitCount(supp)>=2) {
                out[0][3]++; out[1+binID][3]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][4]++; out[1+binID][4]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][5]++; out[1+binID][5]++; }
                if (supp==(MASK_PBSV|MASK_SNIFFLES)) {
                    out[0][9]++; out[1+binID][9]++;
                    if (type==VCFconstants.TYPE_DELETION) { out[0][10]++; out[1+binID][10]++; }
                    else if (type==VCFconstants.TYPE_INSERTION) { out[0][11]++; out[1+binID][11]++; }
                }
                else if (supp==(MASK_PBSV|MASK_PAV)) {
                    out[0][12]++; out[1+binID][12]++;
                    if (type==VCFconstants.TYPE_DELETION) { out[0][13]++; out[1+binID][13]++; }
                    else if (type==VCFconstants.TYPE_INSERTION) { out[0][14]++; out[1+binID][14]++; }
                }
                else if (supp==(MASK_SNIFFLES|MASK_PAV)) {
                    out[0][15]++; out[1+binID][15]++;
                    if (type==VCFconstants.TYPE_DELETION) { out[0][16]++; out[1+binID][16]++; }
                    else if (type==VCFconstants.TYPE_INSERTION) { out[0][17]++; out[1+binID][17]++; }
                }
            }
            if ((supp&MASK_PBSV)!=0) {
                out[0][18]++; out[1+binID][18]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][19]++; out[1+binID][19]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][20]++; out[1+binID][20]++; }
            }
            if ((supp&MASK_SNIFFLES)!=0) {
                out[0][21]++; out[1+binID][21]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][22]++; out[1+binID][22]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][23]++; out[1+binID][23]++; }
            }
            if ((supp&MASK_PAV)!=0) {
                out[0][24]++; out[1+binID][24]++;
                if (type==VCFconstants.TYPE_DELETION) { out[0][25]++; out[1+binID][25]++; }
                else if (type==VCFconstants.TYPE_INSERTION) { out[0][26]++; out[1+binID][26]++; }
            }
            
            // Next iteration
            str=br.readLine();
        }
        br.close();
        
        // Outputting
        System.out.print("-1,"+out[0][0]);
        for (j=1; j<out[0].length; j++) System.out.print(","+out[0][j]);
        System.out.println();
        for (i=0; i<svlengths.length; i++) {
            System.out.print(svlengths[i]+","+out[1+i][0]);
            for (j=1; j<out[1+i].length; j++) System.out.print(","+out[1+i][j]);
            System.out.println();
        }
    }
    
}