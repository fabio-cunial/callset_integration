import java.io.*;


/**
 * Finds windows where hapestry_d32_SVs_and_SNP's haps have the same avg.
 * identity but lower parsimony than hapestry_d3_SVs_and_SNP's haps.
 */
public class CoverageStats6 {

	public static void main(String[] args) throws IOException {
	    final String WINDOW_DIR = args[0];
        
        final String HAPESTRY32_ID = "hapestry_d32_SVs_and_SNPs";
        final String HAPESTRY3_ID = "hapestry_d3_SVs_and_SNPs";
        final String HAPS_FILE = "haps.csv";
        final String NODES_FILE = "nodes.csv";
        
        final double COVERAGE_THRESHOLD = 0.95;
        
        int number;
        int nonrefNodes_hapestry32, nonrefNodes_hapestry3;
        double avgIdentity_hapestry32, avgIdentity_hapestry3, nonrefNodesFullyCovered_hapestry32, nonrefNodesFullyCovered_hapestry3;
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
        
        // Hap identity, 3.
        avgIdentity_hapestry3=0; number=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+HAPESTRY3_ID+"/"+HAPS_FILE)); }
        catch (Exception e) { return; }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            avgIdentity_hapestry3+=Double.parseDouble(tokens[5]);
            number++;
            str=br.readLine();
        }
        br.close();
        avgIdentity_hapestry3/=number;
        
        // Node coverage, 3.
        nonrefNodes_hapestry3=0; nonrefNodesFullyCovered_hapestry3=0;
        try { br = new BufferedReader(new FileReader(WINDOW_DIR+"/"+HAPESTRY3_ID+"/"+NODES_FILE)); }
        catch (Exception e) { 
            System.err.println("File not found: "+WINDOW_DIR+"/"+HAPESTRY3_ID+"/"+NODES_FILE);
            return;
        }
        str=br.readLine();  // Skipping header
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            if (Integer.parseInt(tokens[2])==0) {
                nonrefNodes_hapestry3++;
                if (Double.parseDouble(tokens[4])>=COVERAGE_THRESHOLD) nonrefNodesFullyCovered_hapestry3++;
            }
            str=br.readLine();
        }
        br.close();
        nonrefNodesFullyCovered_hapestry3/=nonrefNodes_hapestry3;
        
        
        
        //System.err.print(WINDOW_DIR+" | "+(avgIdentity_hapestry32-avgIdentity_hapestry3));
        //if (nonrefNodes_hapestry32!=nonrefNodes_hapestry3) System.err.println(WINDOW_DIR+" | "+(nonrefNodesFullyCovered_hapestry32-nonrefNodesFullyCovered_hapestry3)+" | "+nonrefNodes_hapestry32+" | "+nonrefNodes_hapestry3);
        
        //System.err.print(WINDOW_DIR+" | "+avgIdentity_hapestry32+" | "+avgIdentity_hapestry3+" | "+nonrefNodesFullyCovered_hapestry32+" | "+nonrefNodesFullyCovered_hapestry3+" | "+(nonrefNodesFullyCovered_hapestry3/nonrefNodesFullyCovered_hapestry32));
        
        if (avgIdentity_hapestry32>=avgIdentity_hapestry3 && nonrefNodesFullyCovered_hapestry32<nonrefNodesFullyCovered_hapestry3) System.err.println(WINDOW_DIR+" | "+avgIdentity_hapestry32+" | "+avgIdentity_hapestry3+" | "+nonrefNodesFullyCovered_hapestry32+" | "+nonrefNodesFullyCovered_hapestry3+" || "+(avgIdentity_hapestry32/avgIdentity_hapestry3)+" | "+(nonrefNodesFullyCovered_hapestry3/nonrefNodesFullyCovered_hapestry32));
        
        
        //else System.err.println();
        
        //else System.err.println(WINDOW_DIR+" avgIdentity_hapestry32="+avgIdentity_hapestry32+" avgIdentity_hapestry3="+avgIdentity_hapestry3+" avgCoverage_hapestry32="+avgCoverage_hapestry32+" avgCoverage_hapestry3="+avgCoverage_hapestry3);
	}

}