import java.io.*;


/**
 * 
 */
public class BaseQualityDistribution {

	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
        
        final String FASTQ_QUALITIES = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
        
        char c;
        int i, j;
        int nReads;
        String str;
        BufferedReader br;
        long[] histogram;
        
        histogram = new long[FASTQ_QUALITIES.length()];
        br = new BufferedReader(new FileReader(INPUT_FILE));
        nReads=0;
        while (true) {
            str=br.readLine(); 
            if (str==null) break;
            str=br.readLine(); str=br.readLine(); str=br.readLine();
System.err.println(str);            
            for (i=0; i<str.length(); i++) {
                c=str.charAt(i);
                for (j=0; j<FASTQ_QUALITIES.length(); j++) {
                    if (c==FASTQ_QUALITIES.charAt(j)) { histogram[j]++; break; }
                }
            }
            nReads++;
            if (nReads%100000==0) System.err.println(nReads+" reads");
        }
        br.close();
        for (i=0; i<FASTQ_QUALITIES.length(); i++) System.out.println(i+","+histogram[i]);
	}

}