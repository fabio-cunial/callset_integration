import java.io.*;


/**
 * Removes all INFO fields except for few selected ones. Sets FILTER and QUAL to 
 * '.'. Moves the original ID to INFO and assigns a new unique ID.
 */
public class CleanIntersampleVcf {
    
    /**
     * @param args 
     * 0: input file;
     * 1: output file; contains only the '#CHROM' line of the header;
     * >=2: list of all and only the INFO fields that must be preserved; every
     * other INFO field is deleted.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final String OUTPUT_VCF = args[1];
        
        final char COMMENT = '#';
        final String MISSING = ".";
        final String ORIGINAL_ID_FIELD = "TRUVARI_ID";
        final String VCF_HEADER = "#CHROM";
        final char ID_SEPARATOR_OLD = ';';
        final char ID_SEPARATOR_NEW = '_';
        
        int i;
        int nCalls;
        String str, field;
        StringBuilder sb;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        
        sb = new StringBuilder();
        nCalls=0;
        bw = new BufferedWriter(new FileWriter(OUTPUT_VCF));
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) { 
                if (str.substring(0,VCF_HEADER.length()).equals(VCF_HEADER)) { bw.write(str); bw.newLine(); }
                str=br.readLine();
                continue;
            }
            nCalls++;
            if (nCalls%100000==0) System.err.println("Processed "+nCalls+" calls");
            tokens=str.split("\t");
            
            // INFO
            sb.delete(0,sb.length()); sb.append(ORIGINAL_ID_FIELD+"="+tokens[2].replace(ID_SEPARATOR_OLD,ID_SEPARATOR_NEW)+";");
            for (i=2; i<args.length; i++) {
                field=VCFconstants.getField(tokens[7],args[i]);
                if (field!=null) sb.append(args[i]+"="+field+";");
            }
            tokens[7]=sb.toString();
            
            // ID, QUAL, FILTER
            tokens[2]=""+(nCalls-1);
            tokens[5]=MISSING;
            tokens[6]=MISSING;
            
            // Outputting
            bw.write(tokens[0]);
            for (i=1; i<tokens.length; i++) bw.write("\t"+tokens[i]);
            bw.newLine();
            
            // Next iteration
            str=br.readLine();
        }
        bw.close();
    }
    
}