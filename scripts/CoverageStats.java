import java.io.*;


/**
 * Finds windows where truvari's haps have greater avg. coverage than bcftools
 * merge's haps.
 */
public class CoverageStats {

	public static void main(String[] args) throws IOException {
	    final String WINDOW_DIR = args[0];
        
        final String BCFTOOLS_ID = "bcftools_merge_chr1_norm";
        final String TRUVARI_ID = "truvari_chr1_norm";
        final String HAPS_FILE = "haps.csv";
        
        int number;
        double avgBcfTools, avgTruvari;
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
            avgBcfTools+=Double.parseDouble(tokens[4]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgBcfTools/=number;
        
        avgTruvari=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+TRUVARI_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { return; }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgTruvari+=Double.parseDouble(tokens[4]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgTruvari/=number;
        
        if (avgTruvari>avgBcfTools) System.err.println(WINDOW_DIR+" avgTruvari="+avgTruvari+" > avgBcfTools="+avgBcfTools);
	}

}