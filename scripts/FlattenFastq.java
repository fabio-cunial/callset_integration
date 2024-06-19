import java.io.*;


/**
 * Concatenates the four lines of every read.
 *
 * To undo: sed 's/SEPARATOR/\n/g' flattened.txt
 */
public class FlattenFastq {

	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
        final String SEPARATOR = args[1];
        final String OUTPUT_FILE = args[2];
        int BUFFER_SIZE;
        try { BUFFER_SIZE=Integer.parseInt(args[3]); } catch(Exception e) { BUFFER_SIZE=1000000000; }
        
        int nReads;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE),BUFFER_SIZE);
        br = new BufferedReader(new FileReader(INPUT_FILE),BUFFER_SIZE);
        str=br.readLine(); nReads=0;
        while (str!=null) {
            bw.write(str); bw.write(SEPARATOR);
            bw.write(br.readLine()); bw.write(SEPARATOR);
            bw.write(br.readLine()); bw.write(SEPARATOR);
            bw.write(br.readLine()); bw.newLine();
            nReads++;
            if (nReads%100000==0) System.err.println(nReads+" reads");
            str=br.readLine();
        }
        br.close(); bw.close();
        System.err.println("Done");
	}

}