import java.io.*;


/**
 * Finds windows where truvari's graph has zero node coverage.
 */
public class CoverageStats4 {

	public static void main(String[] args) throws IOException {
	    final String WINDOW_DIR = args[0];
        
        final String TRUVARI_ID = "truvari";
        final String HAPS_FILE = "nodes.csv";
        
        double sumTruvari;
        String str;
        BufferedReader br;
        String[] tokens;
        
        sumTruvari=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+TRUVARI_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { return; }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            if (Integer.parseInt(tokens[2])==0) sumTruvari+=Double.parseDouble(tokens[4]);
            str=br.readLine();
        }
        br.close();
        System.err.print(sumTruvari+" ");     
        if (sumTruvari==0) { System.err.println("\n"+WINDOW_DIR+" non-ref nodes have zero coverage\n"); }
	}

}