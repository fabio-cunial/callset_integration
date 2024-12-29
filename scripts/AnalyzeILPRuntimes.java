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
        long value;
        long compressed_d, compressed_n_given_d, compressed_d_plus_n;
        long uncompressed_d, uncompressed_n_given_d, uncompressed_d_plus_n;
        long compression_d, compression_n_given_d, compression_d_plus_n, compression_d_next, compression_n_given_d_next, compression_d_plus_n_next;
        long mandatory, hapsGlobal, hapsLocal, reads, samples, mandatory_next, hapsGlobal_next, hapsLocal_next, reads_next, samples_next;
        long hasLargeWeight1, hasLargeWeight2, easySamples, hasLargeWeight1_next, hasLargeWeight2_next, easySamples_next;
		String str;
		BufferedReader br;
        String[] tokens;
        
        inCompressed=false; inUncompressed=false;
        compressed_d=-1; compressed_n_given_d=-1; compressed_d_plus_n=-1;
        compression_d=-1; compression_n_given_d=-1; compression_d_plus_n=-1;
        compression_d_next=-1; compression_n_given_d_next=-1; compression_d_plus_n_next=-1;
        mandatory=-1; hapsGlobal=-1; hapsLocal=-1; reads=-1; samples=-1;
        mandatory_next=-1; hapsGlobal_next=-1; hapsLocal_next=-1; reads_next=-1; samples_next=-1;
        hasLargeWeight1=-1; hasLargeWeight2=-1; easySamples=-1; hasLargeWeight1_next=-1; hasLargeWeight2_next=-1; easySamples_next=-1;
        uncompressed_d=-1; uncompressed_n_given_d=-1; uncompressed_d_plus_n=-1;
        compressed_n_timeouts=0; uncompressed_n_timeouts=0;
        nEdges_d=0; nEdges_n_given_d=0; nEdges_d_plus_n=0;
        br = new BufferedReader(new FileReader(INPUT_LOG));
		str=br.readLine();
        while (str!=null) {
            System.err.println("Processing line "+str);            
            if (str.indexOf(COMPRESSED_FLAG)==0) {
                // Printing stats about the previous window
                if (compressed_d==-1 || compressed_n_given_d==-1 || compressed_d_plus_n==-1) compressed_n_timeouts++;
                if (uncompressed_d==-1 || uncompressed_n_given_d==-1 || uncompressed_d_plus_n==-1) uncompressed_n_timeouts++;
                if ( compressed_d!=-1 && compressed_n_given_d!=-1 && compressed_d_plus_n!=-1 && 
                     uncompressed_d!=-1 && uncompressed_n_given_d!=-1 && uncompressed_d_plus_n!=-1
                   ) {
                       if (compression_d==-1 || compression_n_given_d==-1 || compression_d_plus_n==-1) {
                           System.err.println("Error: a compression time is -1");
                           System.exit(1);
                       }       
                       System.out.println(compressed_d+","+compression_d+","+compressed_n_given_d+","+compression_n_given_d+","+compressed_d_plus_n+","+compression_d_plus_n+","+    uncompressed_d+","+uncompressed_n_given_d+","+uncompressed_d_plus_n+","+    nEdges_d+","+nEdges_n_given_d+","+nEdges_d_plus_n+","+    mandatory+","+hapsGlobal+","+hapsLocal+","+reads+","+samples+","+(hasLargeWeight1+hasLargeWeight2)+","+easySamples);
                }
                compression_d=compression_d_next; compression_n_given_d=compression_n_given_d_next; compression_d_plus_n=compression_d_plus_n_next;
                mandatory=mandatory_next; hapsGlobal=hapsGlobal_next; hapsLocal=hapsLocal_next; reads=reads_next; samples=samples_next; hasLargeWeight1=hasLargeWeight1_next; hasLargeWeight2=hasLargeWeight2_next; easySamples=easySamples_next;
                compression_d_next=-1; compression_n_given_d_next=-1; compression_d_plus_n_next=-1;
                mandatory_next=-1; hapsGlobal_next=-1; hapsLocal_next=-1; reads_next=-1; samples_next=-1; hasLargeWeight1_next=-1; hasLargeWeight2_next=-1; easySamples_next=-1;
                inCompressed=true; inUncompressed=false;
                compressed_d=-1; compressed_n_given_d=-1; compressed_d_plus_n=-1;
                nEdges_d=0; nEdges_n_given_d=0; nEdges_d_plus_n=0;
            }
            else if (str.indexOf(UNCOMPRESSED_FLAG)==0) {
                inUncompressed=true; inCompressed=false;
                uncompressed_d=-1; uncompressed_n_given_d=-1; uncompressed_d_plus_n=-1;
            }
            else if (str.indexOf(OPTIMIZE_D_PLUS_N_FLAG)==0) {
                tokens=str.split(",");
                value=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                if (inCompressed) compressed_d_plus_n=value;
                else { 
                    uncompressed_d_plus_n=value;
                    nEdges_d_plus_n=Integer.parseInt(tokens[6].substring(tokens[6].lastIndexOf("=")+1));
                }
            }
            else if (str.indexOf(OPTIMIZE_D_FLAG)==0) {
                tokens=str.split(",");
                value=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                if (inCompressed) compressed_d=value;
                else {
                    uncompressed_d=value;
                    nEdges_d=Integer.parseInt(tokens[6].substring(tokens[6].lastIndexOf("=")+1));
                }
            }
            else if (str.indexOf(OPTIMIZE_N_GIVEN_D_FLAG)==0) {
                tokens=str.split(",");
                value=Long.parseLong(tokens[1])*60*60*1000+Integer.parseInt(tokens[2])*60*1000+Integer.parseInt(tokens[3])*1000+Integer.parseInt(tokens[4]);
                if (inCompressed) compressed_n_given_d=value;
                else {
                    uncompressed_n_given_d=value;
                    nEdges_n_given_d=Integer.parseInt(tokens[6].substring(tokens[6].lastIndexOf("=")+1));
                }
            }
            else if (str.indexOf(COMPRESSION_TIME_D)==0) {
                p=str.lastIndexOf("ms");
                compression_d_next=Integer.parseInt(str.substring(COMPRESSION_TIME_D.length()+1,p-1));
            }
            else if (str.indexOf(COMPRESSION_TIME_N_GIVEN_D)==0) {
                p=str.lastIndexOf("ms");
                compression_n_given_d_next=Integer.parseInt(str.substring(COMPRESSION_TIME_N_GIVEN_D.length()+1,p-1));
            }
            else if (str.indexOf(COMPRESSION_TIME_D_PLUS_N)==0) {
                tokens=str.substring(COMPRESSION_TIME_D_PLUS_N.length()+1).split(",");
                mandatory_next=Integer.parseInt(tokens[1]);
                hapsGlobal_next=Integer.parseInt(tokens[2]);
                hapsLocal_next=Integer.parseInt(tokens[3]);
                reads_next=Integer.parseInt(tokens[4]);
                samples_next=Integer.parseInt(tokens[5]);
                hasLargeWeight1_next=Integer.parseInt(tokens[6]);
                hasLargeWeight2_next=Integer.parseInt(tokens[7]);
                easySamples_next=Integer.parseInt(tokens[8]);
                compression_d_plus_n_next=mandatory_next+hapsGlobal_next+hapsLocal_next+reads_next+samples_next+hasLargeWeight1_next+hasLargeWeight2_next+easySamples_next;
            }
            str=br.readLine();
        }
        // Last window
        if ( compressed_d!=-1 && compressed_n_given_d!=-1 && compressed_d_plus_n!=-1 && 
             uncompressed_d!=-1 && uncompressed_n_given_d!=-1 && uncompressed_d_plus_n!=-1
           ) System.out.println(compressed_d+","+compression_d+","+compressed_n_given_d+","+compression_n_given_d+","+compressed_d_plus_n+","+compression_d_plus_n+","+    uncompressed_d+","+uncompressed_n_given_d+","+uncompressed_d_plus_n+","+    nEdges_d+","+nEdges_n_given_d+","+nEdges_d_plus_n+","+    mandatory+","+hapsGlobal+","+hapsLocal+","+reads+","+samples+","+(hasLargeWeight1+hasLargeWeight2)+","+easySamples);
	    br.close();
        System.err.println("Windows with a timeout: uncompressed="+(uncompressed_n_timeouts-1)+" compressed="+(compressed_n_timeouts-1));
	}

}