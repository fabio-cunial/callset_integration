import java.io.*;


/**
 * 
 */
public class ClonesCreateSubsets {

	public static void main(String[] args) throws IOException {
		final String INPUT_CSV = args[0];
        final int N_ROWS = Integer.parseInt(args[1]);  // Total
        final int N_ROWS_PREFIX = Integer.parseInt(args[2]);  // Real
        final int QUANTUM = Integer.parseInt(args[3]);
        final String OUTPUT_CSV = args[4];
        
        int i, j;
        int prefixID, last;
        String str, column1, column2;
        BufferedReader br;
        BufferedWriter bw;
        String[][] rows;
        
        // Loading all rows
        rows = new String[N_ROWS][0];
        br = new BufferedReader(new FileReader(INPUT_CSV));
        str=br.readLine();
        for (i=0; i<N_ROWS; i++) rows[i]=br.readLine().split("\t");
        
        // Creating prefix sets
        bw = new BufferedWriter(new FileWriter(OUTPUT_CSV));
        prefixID=0;
        for (i=QUANTUM-1; i<N_ROWS_PREFIX; i+=QUANTUM) {
            bw.write("prefix_"+prefixID+"\t");
            column1="[\""+rows[0][1]+"\""; 
            column2="[\""+rows[0][2]+"\"";
            last=Math.min(i,N_ROWS_PREFIX-1);
            for (j=1; j<=last; j++) {
                column1+=",\""+rows[j][1]+"\"";
                column2+=",\""+rows[j][2]+"\"";
            }
            bw.write(column1+"]\t"+column2+"]\n");
            prefixID++;
        }
        if (N_ROWS_PREFIX-1>i-QUANTUM) {
            bw.write("prefix_"+prefixID+"\t");
            column1="[\""+rows[0][1]+"\""; 
            column2="[\""+rows[0][2]+"\"";
            for (j=1; j<N_ROWS_PREFIX; j++) {
                column1+=",\""+rows[j][1]+"\"";
                column2+=",\""+rows[j][2]+"\"";
            }
            bw.write(column1+"]\t"+column2+"]\n");
            prefixID++;
        }
        for (i=N_ROWS_PREFIX-1+QUANTUM; i<N_ROWS; i+=QUANTUM) {
            bw.write("prefix_"+prefixID+"\t");
            column1="[\""+rows[0][1]+"\""; 
            column2="[\""+rows[0][2]+"\"";
            last=Math.min(i,N_ROWS-1);
            for (j=1; j<=last; j++) {
                column1+=",\""+rows[j][1]+"\"";
                column2+=",\""+rows[j][2]+"\"";
            }
            bw.write(column1+"]\t"+column2+"]\n");
            prefixID++;
        }
        bw.close();
	}

}