import java.io.*;


/**
 * Finds windows where hapestry_d32_SVs_and_SNP's haps have lower avg. identity 
 * than bcftools's haps.
 */
public class CoverageStats7 {

	public static void main(String[] args) throws IOException {
	    final String WINDOW_DIR = args[0];
        
        final String HAPESTRY32_ID = "hapestry_d32_SVs_and_SNPs";
        final String BCFTOOLS_ID = "SVs_and_SNPs";
        final String HAPS_FILE = "haps.csv";
        final String NODES_FILE = "nodes.csv";
        
        final double COVERAGE_THRESHOLD = 0.95;
        
        int number;
        int nonrefNodes_hapestry32, nonrefNodes_bcftools;
        double avgIdentity_hapestry32, avgIdentity_bcftools, nonrefNodesFullyCovered_hapestry32, nonrefNodesFullyCovered_bcftools;
        String str;
        BufferedReader br;
        String[] tokens;
        
        // Hap identity, 32.
        avgIdentity_hapestry32=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+HAPESTRY32_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { 
            System.err.println("File not found: "+WINDOW_DIR+"/"+HAPESTRY32_ID+"/"+HAPS_FILE);
            return; 
        }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgIdentity_hapestry32+=Double.parseDouble(tokens[5]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgIdentity_hapestry32/=number;
        
        // Node coverage, 32.
        nonrefNodes_hapestry32=0; nonrefNodesFullyCovered_hapestry32=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+HAPESTRY32_ID+"/"+NODES_FILE)); }
        catch (Exception e) { 
            System.err.println("File not found: "+WINDOW_DIR+"/"+HAPESTRY32_ID+"/"+NODES_FILE);
            return;
        }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            if (Integer.parseInt(tokens[2])==0) {
                nonrefNodes_hapestry32++;
                if (Double.parseDouble(tokens[4])>=COVERAGE_THRESHOLD) nonrefNodesFullyCovered_hapestry32++;
            }
            str=br.readLine();
        }
        br.close();
        nonrefNodesFullyCovered_hapestry32/=nonrefNodes_hapestry32;
        
        // Hap identity, bcftools.
        avgIdentity_bcftools=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+BCFTOOLS_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { return; }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgIdentity_bcftools+=Double.parseDouble(tokens[5]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgIdentity_bcftools/=number;
        
        // Node coverage, bcftools.
        nonrefNodes_bcftools=0; nonrefNodesFullyCovered_bcftools=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+BCFTOOLS_ID+"/"+NODES_FILE)); }
        catch (Exception e) { 
            System.err.println("File not found: "+WINDOW_DIR+"/"+BCFTOOLS_ID+"/"+NODES_FILE);
            return;
        }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            if (Integer.parseInt(tokens[2])==0) {
                nonrefNodes_bcftools++;
                if (Double.parseDouble(tokens[4])>=COVERAGE_THRESHOLD) nonrefNodesFullyCovered_bcftools++;
            }
            str=br.readLine();
        }
        br.close();
        nonrefNodesFullyCovered_bcftools/=nonrefNodes_bcftools;
        
        if (avgIdentity_hapestry32<avgIdentity_bcftools) System.out.println(WINDOW_DIR+" | "+avgIdentity_hapestry32+" | "+avgIdentity_bcftools+" | "+nonrefNodesFullyCovered_hapestry32+" | "+nonrefNodesFullyCovered_bcftools+" || "+(avgIdentity_bcftools-avgIdentity_hapestry32)+" | "+(nonrefNodesFullyCovered_bcftools-nonrefNodesFullyCovered_hapestry32));
	}

}