import java.io.*;


public class FixBedBoundaries {

	public static void main(String[] args) throws IOException {
        final String FAI_FILE = args[0];
        final String BED_FILE = args[1];
        
        int i;
        int start, end, lastChromosome;
        String str, chr;
        BufferedReader br;
        int[] chrLengths;
        String[] chrNames, tokens;
        
        // Loading FAI
        chrNames = new String[100];
        chrLengths = new int[100];
        br = new BufferedReader(new FileReader(FAI_FILE));
        str=br.readLine(); lastChromosome=-1;
        while (str!=null) {
            tokens=str.split("\t");
            lastChromosome++;
            chrNames[lastChromosome]=tokens[0];
            chrLengths[lastChromosome]=Integer.parseInt(tokens[1]);
            str=br.readLine();
        }
        br.close();
        
        // Fixing BED
        br = new BufferedReader(new FileReader(BED_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            chr=tokens[0];
            start=Integer.parseInt(tokens[1]);
            end=Integer.parseInt(tokens[2]);
            for (i=0; i<=lastChromosome; i++) {
                if (chrNames[i].equalsIgnoreCase(chr)) {
                    if (start<0) start=0;
                    if (end>chrLengths[i]) end=chrLengths[i];
                    break;
                }
            }
            System.out.println(chr+"\t"+start+"\t"+end);            
            str=br.readLine();
        }
        br.close();
	}

}