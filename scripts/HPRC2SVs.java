import java.io.*;


/**
 * Extracts just the SVs from the joint HPRC VCF that contains both SNPs and 
 * SVs.
 */
public class HPRC2SVs {
    
	public static final char COMMENT = '#';
    public static final String ALT_SEPARATOR = ",";
    
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final int MIN_SV_LENGTH=Integer.parseInt(args[1]);
        final String OUTPUT_VCF = args[2];
        
        int i, p;
        int hap1, hap2, nCalls;
        String str, gt;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens, tokensPrime;
        
        br = new BufferedReader(new FileReader(INPUT_VCF));
        bw = new BufferedWriter(new FileWriter(OUTPUT_VCF));
        str=br.readLine(); nCalls=0;
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            if (nCalls%100000==0) System.err.println("Processed "+nCalls+" calls");
            nCalls++;
            tokens=str.split("\t");
            gt=tokens[9]; 
            if (gt.length()==1) {
                hap1=gt.charAt(0)=='.'?-1:Integer.parseInt(gt);
                hap2=hap1;
            }
            else {
                p=gt.indexOf("|");
                if (p<0) p=gt.indexOf("/");
                if (p<0) {
                    hap1=gt.charAt(0)=='.'?-1:Integer.parseInt(gt);
                    hap2=hap1;
                }
                else {
                    hap1=gt.charAt(0)=='.'?-1:Integer.parseInt(gt.substring(0,p)); 
                    hap2=gt.charAt(p+1)=='.'?-1:Integer.parseInt(gt.substring(p+1));
                }
            }
            if (tokens[4].indexOf(ALT_SEPARATOR)>=0) tokensPrime=tokens[4].split(ALT_SEPARATOR);
            else tokensPrime = new String[] {tokens[4]};
            if ( (hap1>0 && (tokensPrime[hap1-1].length()>=MIN_SV_LENGTH || tokens[3].length()-1>=MIN_SV_LENGTH)) ||
                 (hap2>0 && (tokensPrime[hap2-1].length()>=MIN_SV_LENGTH || tokens[3].length()-1>=MIN_SV_LENGTH))
               ) { bw.write(str); bw.newLine(); }
            str=br.readLine();
        }
        br.close(); bw.close();
    }
    
}