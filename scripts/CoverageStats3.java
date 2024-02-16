import java.io.*;


/**
 * Finds windows where dipcall's haps have lower avg. identity than bcftools
 * merge's haps.
 */
public class CoverageStats3 {

	public static void main(String[] args) throws IOException {
	    final String WINDOW_DIR = args[0];
        
        final String BCFTOOLS_ID = "bcftools_merge_chr1_norm";
        final String DIPCALL_ID = "hprc_47_dipcall_chr1_norm";
        final String HAPS_FILE = "haps.csv";
        
        int number;
        double avgBcfTools, avgDipcall;
        String str;
        BufferedReader br;
        String[] tokens;
        
        avgBcfTools=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+BCFTOOLS_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { return; }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgBcfTools+=Double.parseDouble(tokens[5]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgBcfTools/=number;
        
        avgDipcall=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+DIPCALL_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { return; }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgDipcall+=Double.parseDouble(tokens[5]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgDipcall/=number;
        
        if (avgDipcall<avgBcfTools) System.err.println(WINDOW_DIR+" avgDipcall="+avgDipcall+" < avgBcfTools="+avgBcfTools+" ratio="+(avgBcfTools/avgDipcall));
	}

}