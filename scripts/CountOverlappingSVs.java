import java.util.Arrays;
import java.util.zip.GZIPInputStream;

import java.io.IOException;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.BufferedReader;


/**
 * 
 */
public class CountOverlappingSVs {
    /**
     * Variant types
     */
    private static final int DEL = 0;
    private static final int INV = 1;
    private static final int DUP = 2;
    private static final int INS = 3;
    private static final int REPLACEMENT = 4;
    
    /**
     * A maximal set of overlapping intervals.
     *
     * Remark: the last interval in $window$ might not overlap with the previous
     * intervals, but it might be the first interval of the next set of 
     * overlapping intervals.
     */
    private static Interval[] window;
    
    /**
     * Last element in $window$.
     */
    private static int windowLast;
    
    /**
     * Last position of an overlapping interval in $window$ (one-based,
     * inclusive). Might be smaller than the last position of the last interval
     * in $window$, if such an interval does not overlap with the previous ones.
     */
    private static int windowLastPos;
    

    /**
     * 
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        
        int i;
        int nRecords;
        BufferedReader brVCF;
        long[] histogram;
        
        window = new Interval[100];
        for (i=0; i<window.length; i++) window[i] = new Interval();
        histogram = new long[10000];
        windowLast=-1; windowLastPos=-1;
        brVCF = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        nRecords=0;
        while (true) {
            loadNextWindow(brVCF);
            if (windowLast==-1) break;
            countWindow(histogram);
            nRecords+=(isLastInWindow()?windowLast:windowLast-1)+1;
            if (nRecords%10000==0) System.out.println("Loaded "+nRecords+" records");
        }
        brVCF.close();
        for (i=0; i<histogram.length; i++) System.out.println(i+","+histogram[i]);
    }
    
    
    /**
     * Loads in $window$ a maximal set of overlapping intervals, regardless of
     * their genotypes.
     *
     * Remark: after this procedure completes, the last interval in $window$
     * might not overlap with the previous ones.
     */
    private static final void loadNextWindow(BufferedReader br) throws IOException {
        int i;
        String str;
        Interval tmpInterval;
        
        if (isLastInWindow()) {
            windowLast=-1;
            windowLastPos=-1;
        }
        else {
            tmpInterval=window[0];
            window[0]=window[windowLast];
            window[windowLast]=tmpInterval;
            windowLastPos=window[0].last;
            windowLast=0;
        }
        while (true) {
            str=br.readLine();
            if (str==null) break;
            if (str.charAt(0)==COMMENT) continue;
            windowLast++;
            if (windowLast==window.length) {
                Interval[] newArray = new Interval[window.length<<1];
                System.arraycopy(window,0,newArray,0,window.length);
                for (i=window.length; i<newArray.length; i++) newArray[i] = new Interval();
                window=newArray;
            }
            tmpInterval=window[windowLast];
            tmpInterval.init(str);
            if (!isLastInWindow()) break;
            if (tmpInterval.last>windowLastPos) windowLastPos=tmpInterval.last;
        }
    }
    
    
    /**
     * @return TRUE iff the last interval of $window$ overlaps with the previous
     * intervals.
     */
    private static final boolean isLastInWindow() {
        return windowLast==-1 || windowLastPos==-1 || (window[windowLast].chr==window[0].chr && window[windowLast].first<=windowLastPos);
    }
    
    
    /**
     * @param histogram for each $i$, the number of records with $i$ overlaps.
     */
    private static final void countWindow(long[] histogram) throws IOException {
        final int LAST = isLastInWindow()?windowLast:windowLast-1;
        
        int i, j;
        int firstPos, lastPos, nOverlaps;
        
        Interval.order=Interval.ORDER_FIRST_POS;
        Arrays.sort(window,0,LAST+1);
        for (i=0; i<=LAST; i++) {
            firstPos=window[i].first; lastPos=window[i].last;
            nOverlaps=0;
            for (j=i+1; j<=LAST; j++) {
                if (window[j].first>lastPos) break;
                nOverlaps++;
            }
            for (j=i-1; j>=0; j--) {
                if (window[j].last>firstPos) nOverlaps++;
            }
            histogram[nOverlaps>=histogram.length?histogram.length-1:nOverlaps]++;
        }
    }
    
        
    
    
    // -------------------------- DATA STRUCTURES  -----------------------------
    
    /**
     * A VCF record represented as a weighted interval on the line.
     * The same object is intended to be reused by multiple VCF records.
     */
    private static class Interval implements Comparable {
        public static final int ORDER_FIRST_POS = 0;
        public static final int INS_DELTA = 10;  // Arbitrary
        public static int order;
        
        public int variantType;
        public int chr;
        public int first, last;  // One-based, inclusive.
        public String[] vcfRecord;  // One cell per VCF column
        
        
        public Interval() {
            variantType=-1; chr=-1; first=-1; last=-1;
            vcfRecord=null;
        }
        
        
        /**
         * @param record a VCF record.
         */
        public final void init(String record) {
            int i;
            int pos, length;
            String tmpString;
            
            this.vcfRecord=record.split(VCF_SEPARATOR);
            tmpString=getInfoField(this.vcfRecord[7],SVTYPE_STR);
            if (tmpString!=null && tmpString.length()!=0) variantType=svType2Row(tmpString);
            else variantType=refAlt2Row(this.vcfRecord[3],this.vcfRecord[4]);
            if (variantType==-1) {
                System.err.println("ERROR: this record encodes an unknown variant type: "+record);
                System.exit(1);
            }
            chr=string2contig(this.vcfRecord[0]);
            pos=Integer.parseInt(this.vcfRecord[1]);
            tmpString=getInfoField(this.vcfRecord[7],SVLEN_STR);
            if (tmpString!=null) {
                length=Integer.parseInt(tmpString);
                if (length<0) length=-length;
            }
            else if (variantType==REPLACEMENT) length=this.vcfRecord[3].length()-1;
            else length=Math.max(this.vcfRecord[3].length(),this.vcfRecord[4].length())-1;
            first=-1; last=-1;
            if (variantType==DEL || variantType==INV || variantType==DUP || variantType==REPLACEMENT) { first=pos+1; last=pos+length; }
            else if (variantType==INS) { first=pos+1-INS_DELTA; last=pos+INS_DELTA; }
        }

        
        public String toString() {
            return variantType+", chr"+chr+"["+first+".."+last+"]";
        }
        
        
        public boolean equals(Object other) {
            final Interval otherInterval = (Interval)other;
            return first==otherInterval.first;
        }
        
        
        public int compareTo(Object other) {
            final Interval otherInterval = (Interval)other;
            if (first<otherInterval.first) return -1;
            else if (first>otherInterval.first) return 1;
            else return 0;
        }
        
    }
    
    
	private static final int svType2Row(String type) {
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return DEL;
		else if (type.equalsIgnoreCase(INV_STR)) return INV;
        else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR) ||
                  type.equalsIgnoreCase(CNV_STR)
			    ) return DUP;
        else if ( type.equalsIgnoreCase(INS_STR) ||
                  type.equalsIgnoreCase(INS_ME_STR) ||
                  type.equalsIgnoreCase(INS_NOVEL_STR)
                ) return INS;
		else return -1;
	}


	private static final int refAlt2Row(String ref, String alt) {
        if (ref.length()==1) {
            if (alt.length()>1) return INS;
            else return -1;
        }
        else {
            if (alt.length()==1) return DEL;
            else return REPLACEMENT;
        }
	}
    
    
    
    
    // ------------------------- BASIC VCF HANDLING ----------------------------
    
    private static final String VCF_SEPARATOR = "\t";
    private static final char GT_SEPARATOR = ':';
    private static final String GT_STR = "GT";
    public static final char COMMENT = '#';
    public static final String END_STR = "END";
    public static final String CIEND_STR = "CIEND";
    public static final String INFO_SEPARATOR = ";";
    public static final String SVTYPE_STR = "SVTYPE";
    public static final String SVLEN_STR = "SVLEN";
    public static final String CHR_STR = "chr";
    public static final int CHR_STR_LENGTH = CHR_STR.length();
    public static final String X_STR_PRIME = "X";
    public static final String Y_STR_PRIME = "Y";
    public static final String M_STR_PRIME = "M";
    public static final String MT_STR_PRIME = "MT";
    public static final String X_STR = CHR_STR+X_STR_PRIME;
    public static final String Y_STR = CHR_STR+Y_STR_PRIME;
    public static final String M_STR = CHR_STR+M_STR_PRIME;
    public static final String MT_STR = CHR_STR+MT_STR_PRIME;
    
    /**
     * SV types: labels used by callers.
     */
    public static final String DEL_STR = "DEL";
    public static final String DEL_ME_STR = "DEL:ME";
    public static final String DEL_INV_STR = "DEL/INV";
    public static final String INS_STR = "INS";
    public static final String INS_ME_STR = "INS:ME";
    public static final String INS_NOVEL_STR = "INS:NOVEL";
    public static final String DUP_STR = "DUP";
    public static final String DUP_TANDEM_STR = "DUP:TANDEM";
    public static final String DUP_INT_STR = "DUP:INT";
    public static final String INV_STR = "INV";
    public static final String INV_DUP_STR = "INVDUP";
    public static final String CNV_STR = "CNV";
    public static final String BND_STR = "BND";
    public static final String TRA_STR = "TRA";
    
    
    /**
     * @return NULL if $field$ does not occur in $str$.
     */
    public static final String getInfoField(String str, String field) {
        final int FIELD_LENGTH = field.length()+1;
        int p = str.indexOf(field+"=");
        if (p<0) return null;
        if (field.equalsIgnoreCase(END_STR)) {
            while (p>=2 && str.substring(p-2,p-2+CIEND_STR.length()).equalsIgnoreCase(CIEND_STR)) p=str.indexOf(field+"=",p+1);
            if (p<0) return null;
        }
        final int q = str.indexOf(INFO_SEPARATOR,p+FIELD_LENGTH);
        return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
    }
    
    
    /**
     * @return one-based.
     */
    public static final int string2contig(String str) {
        if (str.length()>=CHR_STR_LENGTH && str.substring(0,CHR_STR_LENGTH).equalsIgnoreCase(CHR_STR)) {
            if (str.equalsIgnoreCase(X_STR)) return 23;
            else if (str.equalsIgnoreCase(Y_STR)) return 24;
            else if (str.equalsIgnoreCase(M_STR) || str.equalsIgnoreCase(MT_STR)) return 25;
            else if (str.length()<=CHR_STR_LENGTH+2) return Integer.parseInt(str.substring(CHR_STR_LENGTH));
            else return -1;
        }
        else {
            if (str.equalsIgnoreCase(X_STR_PRIME)) return 23;
            else if (str.equalsIgnoreCase(Y_STR_PRIME)) return 24;
            else if (str.equalsIgnoreCase(M_STR_PRIME) || str.equalsIgnoreCase(MT_STR_PRIME)) return 25;
            else if (str.length()<=2) return Integer.parseInt(str);
            else return -1;
        }
    }
    

}