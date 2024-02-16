import java.io.*;


/**
 * Finds windows where jasmine's haps have nonzero avg. coverage, but bcftools
 * merge haps have zero coverage.
 */
public class CoverageStats2 {

	public static void main(String[] args) throws IOException {
	    final String WINDOW_DIR = args[0];
        
        final String BCFTOOLS_ID = "bcftools_merge_chr1_norm";
        final String JASMINE_ID = "jasmine_default_chr1_norm";
        final String HAPS_FILE = "haps.csv";
        
        int number;
        double avgBcfTools, avgJasmine;
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
        
        avgJasmine=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+JASMINE_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { return; }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgJasmine+=Double.parseDouble(tokens[4]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgJasmine/=number;
        
        if (avgBcfTools==0/* && avgJasmine>0*/) System.err.println(WINDOW_DIR+" avgJasmine="+avgJasmine+" > avgBcfTools="+avgBcfTools);
	}

}