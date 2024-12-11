import java.io.*;


public class AnalyzeILPRuntimes {

    /**
     * @param args
     */
	public static void main(String[] args) throws IOException {
		final String INPUT_LOG = args[0];
		
        final String COMPRESSED_FLAG = "---- COMPRESSED";
        final String UNCOMPRESSED_FLAG = "---- UNCOMPRESSED";
        final String OPTIMIZE_D_FLAG = "optimize_d";
        final String OPTIMIZE_N_GIVEN_D_FLAG = "optimize_n_given_d";
        final String OPTIMIZE_D_PLUS_N_FLAG = "optimize_d_plus_n";
        final String COMPRESSION_TIME_D = "Compression time d:";
        final String COMPRESSION_TIME_N_GIVEN_D= "Compression time n_given_d:";
        final String COMPRESSION_TIME_D_PLUS_N = "Compression time d_plus_n:";
        
        boolean inCompressed, inUncompressed;
        int i, p;
        int compressed_n_timeouts, uncompressed_n_timeouts, nEdges_d, nEdges_n_given_d, nEdges_d_plus_n;
        long compressed_d, compressed_n_given_d, compressed_d_plus_n;
        long uncompressed_d, uncompressed_n_given_d, uncompressed_d_plus_n;
        long compression_d, compression_n_given_d, compression_d_plus_n;
        long mandatory, hapsGlobal, hapsLocal, reads, samples;
		String str;
		BufferedReader br;
        String[] tokens;
        
        inCompressed=false; inUncompressed=false;
        compressed_d=-1; compression_d=-1; compressed_n_given_d=-1; compression_n_given_d=-1; compressed_d_plus_n=-1; compression_d_plus_n=-1;
        mandatory=-1; hapsGlobal=-1; hapsLocal=-1; reads=-1; samples=-1;
        uncompressed_d=-1; uncompressed_n_given_d=-1; uncompressed_d_plus_n=-1;
        compressed_n_timeouts=0; uncompressed_n_timeouts=0;
        nEdges_d=0; nEdges_n_given_d=0; nEdges_d_plus_n=0;
        br = new BufferedReader(new FileReader(INPUT_LOG));
		str=br.readLine();
        while (str!=null) {
            if (str.indexOf(COMPRESSED_FLAG)==0) {
                if (compressed_d==-1 || compressed_n_given_d==-1 || compressed_d_plus_n==-1) compressed_n_timeouts++;
                if (uncompressed_d==-1 || uncompressed_n_given_d==-1 || uncompressed_d_plus_n==-1) uncompressed_n_timeouts++;
                if ( compressed_d!=-1 && compressed_n_given_d!=-1 && compressed_d_plus_n!=-1 && 
                     uncompressed_d!=-1 && uncompressed_n_given_d!=-1 && uncompressed_d_plus_n!=-1
                   ) System.out.println(compressed_d+","+compression_d+","+compressed_n_given_d+","+compression_n_given_d+","+compressed_d_plus_n+","+compression_d_plus_n+","+    uncompressed_d+","+uncompressed_n_given_d+","+uncompressed_d_plus_n+","+    nEdges_d+","+nEdges_n_given_d+","+nEdges_d_plus_n+","+    mandatory+","+hapsGlobal+","+hapsLocal+","+reads+","+samples);
                inCompressed=true; inUncompressed=false;
                compressed_d=-1; compression_d=-1; compressed_n_given_d=-1; compression_n_given_d=-1; compressed_d_plus_n=-1; compression_d_plus_n=-1;
                mandatory=-1; hapsGlobal=-1; hapsLocal=-1; reads=-1; samples=-1;
                nEdges_d=0; nEdges_n_given_d=0; nEdges_d_plus_n=0;
            }
            else if (str.indexOf(UNCOMPRESSED_FLAG)==0) {
                inUncompressed=true; inCompressed=false;
                uncompressed_d=-1; uncompressed_n_given_d=-1; uncompressed_d_plus_n=-1;
            }
            else if (str.indexOf(OPTIMIZE_D_PLUS_N_FLAG)==0) {
                tokens=str.split(",");
                if (inCompressed) compressed_d_plus_n=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                else { 
                    uncompressed_d_plus_n=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                    nEdges_d_plus_n=Integer.parseInt(tokens[6].substring(tokens[6].lastIndexOf("=")+1));
                }
            }
            else if (str.indexOf(OPTIMIZE_D_FLAG)==0) {
                tokens=str.split(",");
                if (inCompressed) compressed_d=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                else {
                    uncompressed_d=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                    nEdges_d=Integer.parseInt(tokens[6].substring(tokens[6].lastIndexOf("=")+1));
                }
            }
            else if (str.indexOf(OPTIMIZE_N_GIVEN_D_FLAG)==0) {
                tokens=str.split(",");
                if (inCompressed) compressed_n_given_d=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                else {
                    uncompressed_n_given_d=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                    nEdges_n_given_d=Integer.parseInt(tokens[6].substring(tokens[6].lastIndexOf("=")+1));
                }
            }
            else if (str.indexOf(COMPRESSION_TIME_D)==0) {
                p=str.lastIndexOf("ms");
                compression_d=Integer.parseInt(str.substring(COMPRESSION_TIME_D.length()+1,p-1));
            }
            else if (str.indexOf(COMPRESSION_TIME_N_GIVEN_D)==0) {
                p=str.lastIndexOf("ms");
                compression_n_given_d=Integer.parseInt(str.substring(COMPRESSION_TIME_N_GIVEN_D.length()+1,p-1));
            }
            else if (str.indexOf(COMPRESSION_TIME_D_PLUS_N)==0) {
                tokens=str.substring(COMPRESSION_TIME_D_PLUS_N.length()+1).split(",");
                mandatory=Integer.parseInt(tokens[1]);
                hapsGlobal=Integer.parseInt(tokens[2]);
                hapsLocal=Integer.parseInt(tokens[3]);
                reads=Integer.parseInt(tokens[4]);
                samples=Integer.parseInt(tokens[5]);
                compression_d_plus_n=mandatory+hapsGlobal+hapsLocal+reads+samples;
            }
            str=br.readLine();
        }
	    br.close();
        System.err.println("Windows with a timeout: uncompressed="+(uncompressed_n_timeouts-1)+" compressed="+(compressed_n_timeouts-1));
	}

}