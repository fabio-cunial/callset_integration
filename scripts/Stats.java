import java.io.*;


public class Stats {

	public static void main(String[] args) throws IOException {
	    final String INPUT_FILE = args[0];
        
        final int MIN_SV_LENGTH = 0;
        final int MAX_SV_LENGTH = 5000;
        final int SV_LENGTH_QUANTUM = 10;
        final double FRACTION_QUANTUM = 0.01;
        
        int i, j;
        int length;
        double fraction, sum;
        String str;
        BufferedReader br;
        String[] tokens;
        double[][] matrix;
        
        // Collecting counts
        matrix = new double[1+MAX_SV_LENGTH/SV_LENGTH_QUANTUM][1+(int)(1.0/FRACTION_QUANTUM)];
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            length=Integer.parseInt(tokens[0]);
            if (length>MAX_SV_LENGTH) {
                str=br.readLine();
                continue;
            }
            fraction=Double.parseDouble(tokens[1]);
            matrix[length/SV_LENGTH_QUANTUM][(int)(fraction/FRACTION_QUANTUM)]+=1.0;
            str=br.readLine();
        }
        br.close();
        
        // Normalizing counts
        // for (i=0; i<matrix.length; i++) {
        //     sum=0;
        //     for (j=0; j<matrix[i].length; j++) sum+=matrix[i][j];
        //     if (sum!=0) {
        //         for (j=0; j<matrix[i].length; j++) matrix[i][j]/=sum;
        //     }
        // }
        
        // Printing
        for (i=0; i<matrix.length; i++) {
            System.out.print((i*SV_LENGTH_QUANTUM)+"");
            for (j=0; j<matrix[i].length; j++) System.out.print(","+matrix[i][j]);
            System.out.println();
        }
	}

}