import java.io.*;


/**
 * Finds windows where hapestry_d32_SVs_and_SNP's haps have lower avg. identity 
 * than SVs_and_SNPs's haps.
 */
public class CoverageStats5 {

	public static void main(String[] args) throws IOException {
	    final String WINDOW_DIR = args[0];
        
        final String HAPESTRY_ID = "hapestry_d32_SVs_and_SNPs";
        final String INPUT_ID = "SVs_and_SNPs";
        final String HAPS_FILE = "haps.csv";
        
        int number;
        double avgHapestry, avgInput;
        String str;
        BufferedReader br;
        String[] tokens;
        
        avgHapestry=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+HAPESTRY_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { 
            System.err.println("File not found: "+WINDOW_DIR+"/"+HAPESTRY_ID+"/"+HAPS_FILE);
            return; 
        }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgHapestry+=Double.parseDouble(tokens[5]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgHapestry/=number;
        
        avgInput=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+INPUT_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { return; }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgInput+=Double.parseDouble(tokens[5]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgInput/=number;
        
        if (avgHapestry<avgInput) System.err.println(WINDOW_DIR+" avgHapestry="+avgHapestry+" < avgInput="+avgInput+" ratio="+(avgInput/avgHapestry));
	}

}