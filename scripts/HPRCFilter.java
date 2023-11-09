import java.io.*;


/**
 * Given an HPRC joint VCF and a list of samples, the program creates a 
 * separate truth VCF for each sample.
 *
 * java -Xmx16G  HGSVC2Filter hprc-v1.0-mc-grch38.vcf . 1   2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45
 */
public class HPRCFilter {
    /**
     * The program writes to the VCF of sample X all and only the calls whose GT 
     * in sample X is either *|Y or Y|* with Y>0.
     *
     * @param args 
     * 3: list of one-based, comma-separated column indexes to select, where 
     *    1=first sample column.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final String OUT_DIR = args[1];
        final String SELECTED_SAMPLES = args[2];
        
        int i, j, p;
        int nColumns, hap1, hap2;
        String str, gt;
        StringBuilder header;
        BufferedReader br;
        int[] columns;
        String[] tokens, tokensPrime;
        BufferedWriter[] outputs;
        
        tokens=SELECTED_SAMPLES.split(",");
        nColumns=tokens.length;
        columns = new int[nColumns];
        for (i=0; i<nColumns; i++) columns[i]=Integer.parseInt(tokens[i]);
        outputs = new BufferedWriter[nColumns];
        header = new StringBuilder();
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.substring(0,6).equals("#CHROM")) {
                tokens=str.split("\t");
                for (i=0; i<nColumns; i++) {
                    outputs[i] = new BufferedWriter(new FileWriter(OUT_DIR+"/"+tokens[8+columns[i]]+".vcf"));
                    outputs[i].write(header.toString());
                    for (j=0; j<=8; j++) outputs[i].write(tokens[j]+"\t");
                    outputs[i].write(tokens[8+columns[i]]);
                    outputs[i].newLine();
                }
                str=br.readLine();
                continue;
            }
            else if (str.charAt(0)=='#') {
                header.append(str+"\n");
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            for (i=0; i<nColumns; i++) {
                gt=tokens[8+columns[i]];
                p=gt.indexOf("|");
                if (p<0) p=gt.indexOf("/");
                if (p<0) {
                    hap1=gt.equals(".")?-1:Integer.parseInt(gt);
                    hap2=hap1;
                }
                else {
                    hap1=gt.substring(0,p).equals(".")?-1:Integer.parseInt(gt.substring(0,p));
                    hap2=gt.substring(p+1).equals(".")?-1:Integer.parseInt(gt.substring(p+1));
                }
                if (hap1>0 || hap2>0) {
                    for (j=0; j<=8; j++) outputs[i].write(tokens[j]+"\t");
                    outputs[i].write(tokens[8+columns[i]]);
                    outputs[i].newLine();
                }
            }
            str=br.readLine();
        }
        br.close();
        for (i=0; i<nColumns; i++) outputs[i].close();
    }
    
}