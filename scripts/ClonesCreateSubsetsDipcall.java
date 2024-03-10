import java.io.*;


/**
 * 
 */
public class ClonesCreateSubsetsDipcall {

	public static void main(String[] args) throws IOException {
		final String INPUT_SAMPLES_LIST = args[0];
        final int N_SAMPLES = Integer.parseInt(args[1]);
        final String DIPCALL_REMOTE_DIR = args[2];
        final int QUANTUM = Integer.parseInt(args[3]);
        final String OUTPUT_CSV = args[4];
        
        int i, j;
        int prefixID, last;
        String str, column1, column2;
        BufferedReader br;
        BufferedWriter bw;
        String[] samples;
        
        // Loading samples
        samples = new String[N_SAMPLES];
        br = new BufferedReader(new FileReader(INPUT_SAMPLES_LIST));
        for (i=0; i<N_SAMPLES; i++) samples[i]=br.readLine();
        br.close();
        
        // Creating prefix sets
        bw = new BufferedWriter(new FileWriter(OUTPUT_CSV));
        prefixID=0;
        for (i=QUANTUM-1; i<N_SAMPLES; i+=QUANTUM) {
            bw.write("prefix_"+prefixID+"\t");
            column1="[\""+DIPCALL_REMOTE_DIR+"/"+samples[0]+"_sv.vcf.gz\""; 
            column2="[\""+DIPCALL_REMOTE_DIR+"/"+samples[0]+"_sv.vcf.gz.tbi\"";
            last=Math.min(i,N_SAMPLES-1);
            for (j=1; j<=last; j++) {
                column1+=",\""+DIPCALL_REMOTE_DIR+"/"+samples[j]+"_sv.vcf.gz\""; 
                column2+=",\""+DIPCALL_REMOTE_DIR+"/"+samples[j]+"_sv.vcf.gz.tbi\"";
            }
            bw.write(column1+"]\t"+column2+"]\n");
            prefixID++;
        }
        if (N_SAMPLES-1>i-QUANTUM) {
            bw.write("prefix_"+prefixID+"\t");
            column1="[\""+DIPCALL_REMOTE_DIR+"/"+samples[0]+"_sv.vcf.gz\""; 
            column2="[\""+DIPCALL_REMOTE_DIR+"/"+samples[0]+"_sv.vcf.gz.tbi\"";
            for (j=1; j<N_SAMPLES; j++) {
                column1+=",\""+DIPCALL_REMOTE_DIR+"/"+samples[j]+"_sv.vcf.gz\""; 
                column2+=",\""+DIPCALL_REMOTE_DIR+"/"+samples[j]+"_sv.vcf.gz.tbi\"";
            }
            bw.write(column1+"]\t"+column2+"]\n");
            prefixID++;
        }
        bw.close();
	}

}