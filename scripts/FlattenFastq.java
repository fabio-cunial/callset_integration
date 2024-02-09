import java.io.*;


/**
 * Concatenates the four lines of every read.
 * To undo: sed 's/SEPARATOR/\n/g' flattened.txt
 */
public class FlattenFastq {

	public static void main(String[] args) throws IOException {
		final String INPUT_FILE = args[0];
        final String SEPARATOR = args[1];
        final String OUTPUT_FILE = args[2];
        
        int nReads;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine(); nReads=0;
        while (str!=null) {
            bw.write(str); bw.write(SEPARATOR);
            bw.write(br.readLine()); bw.write(SEPARATOR);
            bw.write(br.readLine()); bw.write(SEPARATOR);
            bw.write(br.readLine()); bw.write(SEPARATOR);
            bw.newLine();
            nReads++;
            if (nReads%100000==0) System.err.println(nReads+" reads");
            str=br.readLine();
        }
        br.close(); bw.close();
        System.err.println("Done");
	}

}