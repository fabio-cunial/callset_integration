import java.io.*;


/**
 * Prepares for sv-merger a TRF file downloaded from e.g. 
 * https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
 */
public class TRF2SVMerger {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_BED = args[0];
        
        int id;
        String str, currentChromosome;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(INPUT_BED));
        bw=null;
        currentChromosome=""; id=0;
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            if (!tokens[0].equalsIgnoreCase(currentChromosome)) {
                if (bw!=null) bw.close();
                currentChromosome=tokens[0];
                bw = new BufferedWriter(new FileWriter(currentChromosome+".trf.sorted.gor"));
                id=0;
            }
            bw.write(tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+(Integer.parseInt(tokens[2])-Integer.parseInt(tokens[1]))+"\t"+(id++));
            bw.newLine();
            str=br.readLine();
        }
        br.close();
    }
    
}