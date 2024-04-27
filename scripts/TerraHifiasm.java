import java.io.*;


/**
 * Reformats a Terra table so that it contains the two assembled haplotypes from
 * hifiasm, rather than their containing directory.
 */
public class TerraHifiasm {
    /**
     * @param args
     * 1: typically "hifiasm_haplotigs";
     * 3,4: new columns that contain the assembled haplotypes from the haplotigs
     * column.
     */
	public static void main(String[] args) throws IOException {
	    final String INPUT_TSV = args[0];
        final String SAMPLE_ID_COLUMN = args[1];
        final String HAPLOTIGS_COLUMN = args[2];
        final String HAP1_COLUMN = args[3];
        final String HAP2_COLUMN = args[4];
        
        final String SEPARATOR = "\t";
        
        int i;
        int haplotigsColumn, idColumn;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(INPUT_TSV));
        str=br.readLine();
        tokens=str.split(SEPARATOR);
        idColumn=-1; haplotigsColumn=-1;
        for (i=0; i<tokens.length; i++) {
            if (tokens[i].equalsIgnoreCase(SAMPLE_ID_COLUMN)) idColumn=i;
            if (tokens[i].equalsIgnoreCase(HAPLOTIGS_COLUMN)) haplotigsColumn=i;
        }
        if (haplotigsColumn==-1) { System.err.println("ERROR: no column "+HAPLOTIGS_COLUMN+" in file "+INPUT_TSV); System.exit(1); }
        if (idColumn==-1) { System.err.println("ERROR: no column "+SAMPLE_ID_COLUMN+" in file "+INPUT_TSV); System.exit(1); }
        System.err.println("HAPLOTIGS_COLUMN="+HAPLOTIGS_COLUMN+"="+haplotigsColumn+" SAMPLE_ID_COLUMN="+SAMPLE_ID_COLUMN+"="+idColumn);
        System.out.println(str+SEPARATOR+HAP1_COLUMN+SEPARATOR+HAP2_COLUMN);  //
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            System.out.println(str+SEPARATOR+tokens[haplotigsColumn]+tokens[idColumn]+".bp.hap1.p_ctg.fa.gz"+SEPARATOR+tokens[haplotigsColumn]+tokens[idColumn]+".bp.hap2.p_ctg.fa.gz");
            str=br.readLine();
        }
        br.close();
	}

}