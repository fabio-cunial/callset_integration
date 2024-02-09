import java.io.*;


/**
 * Given a merged VCF, the program keeps just the calls that were not created 
 * exclusively by PBSV or Sniffles.
 */
public class FilterMergedVcf {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final boolean REMOVE_PBSV = Integer.parseInt(args[1])==1;
        final boolean REMOVE_SNIFFLES = Integer.parseInt(args[2])==1;
        final String OUTPUT_VCF = args[3];
        
        final char COMMENT = '#';
        final String ID_SEPARATOR = ",";
        final String ID_PBSV = "pbsv";
        final String ID_SNIFFLES = "sniffles";
        
        boolean isMerged;
        int nCalls, nMergedCalls, nPbsvCalls, nSnifflesCalls, nCallsKept;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        
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
            tokens[2]=tokens[2].toLowerCase();
            nMergedCalls=nOccurrences(tokens[2],ID_SEPARATOR)+1;
            nPbsvCalls=nOccurrences(tokens[2],ID_PBSV);
            nSnifflesCalls=nOccurrences(tokens[2],ID_SNIFFLES);
            if (nMergedCalls>nPbsvCalls+nSnifflesCalls) { bw.write(str); bw.newLine(); nCallsKept++; }
            str=br.readLine();
        }
        br.close(); bw.close();
        System.err.println("Kept "+nCallsKept+" calls out of "+nCalls+" total ("+()+"%)");
    }

    
    private static final int nOccurrences(String text, String query) {
        int i = 0; 
        int n = 0;
        do {
            i=text.indexOf(query,i);
            if (i<0) break;
            n++; i++;
        } while (i<text.length());
        return n;
    }
}