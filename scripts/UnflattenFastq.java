import java.util.Arrays;
import java.io.*;


/**
 * Undoes $FlattenFastq.java$, and prints to STDOUT the last line (one-based,
 * one value per row) in the unflattened file that corresponds to each coverage
 * provided in input. This can be used with $head$ to create the coverage
 * subsamples.
 */
public class UnflattenFastq {

	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
        final String SEPARATOR = args[1];
        final long GENOME_LENGTH = Long.parseLong(args[2]);
        final String COVERAGES = args[3];  // Integer, comma-separated.
        final String OUTPUT_FILE = args[4];
        int BUFFER_SIZE;
        try { BUFFER_SIZE=Integer.parseInt(args[5]); } catch(Exception e) { BUFFER_SIZE=1000000000; }
        
        final int SEPARATOR_LENGTH = SEPARATOR.length();
        
        int i, j, p, q;
        int nCoverages, nReads;
        long nChars;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        int[] lastLine;
        long[] coverageChars;        
        int[] coverages;
        String[] tokens;
        
        // Computing n. characters per coverage
        tokens=COVERAGES.split(",");
        nCoverages=tokens.length;
        coverages = new int[nCoverages];
        for (i=0; i<nCoverages; i++) coverages[i]=Integer.parseInt(tokens[i]);
        Arrays.sort(coverages);
        coverageChars = new long[nCoverages];
        for (i=0; i<nCoverages; i++) coverageChars[i]=(long)(coverages[i]*GENOME_LENGTH);
        lastLine = new int[nCoverages];
        Arrays.fill(lastLine,-1);
        
        // Unflattening
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE),BUFFER_SIZE);
        br = new BufferedReader(new FileReader(INPUT_FILE),BUFFER_SIZE);
        str=br.readLine(); nReads=0; nChars=0; j=0;
        while (str!=null) {
            nReads++;
            p=str.indexOf(SEPARATOR); bw.write(str,0,p); bw.newLine();
            p+=SEPARATOR_LENGTH; q=str.indexOf(SEPARATOR,p); bw.write(str,p,q-p); bw.newLine();
            nChars+=q-p;
            if (nChars>=coverageChars[j]) lastLine[j++]=nReads*4;  // One-based
            p=q+SEPARATOR_LENGTH; q=str.indexOf(SEPARATOR,p); bw.write(str,p,q-p); bw.newLine();
            p=q+SEPARATOR_LENGTH; bw.write(str,p,str.length()-p); bw.newLine();
            if (nReads%100000==0) System.err.println(nReads+" reads");
            if (j==nCoverages) break;
            str=br.readLine();
        }
        br.close(); bw.close();
        System.err.println("Done");
        
        // Outputting the last line of every coverage (one-based).
        for (i=0; i<nCoverages; i++) System.out.println(coverages[i]+","+(lastLine[i]!=-1?lastLine[i]:(nReads*4)));
        System.out.println();
	}

}