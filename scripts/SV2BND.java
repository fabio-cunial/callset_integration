import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class SV2BND {
    /**
     * Chromosome buffer
     */
    private static String CHROMOSOMES_DIR;
    private static String currentChromosome;  // ID
    private static StringBuilder sb;  // Sequence
    
    
    /**
     * @param args
     * 1: only calls this length or longer are used;
     * 2: converts and INS into a single breakend with this number of bases.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        CHROMOSOMES_DIR=args[1];
        final int MIN_SV_LENGTH = Integer.parseInt(args[2]);
        final int SINGLE_BREAKEND_LENGTH = Integer.parseInt(args[3]);
        final String OUTPUT_VCF = args[4];
        
        final char COMMENT = '#';
        final String MATEID_STR = "MATEID";
        final String NEW_MATE_PREFIX = "mate_of_";
        final String SECOND_ADJACENCY_SUFFIX = "_prime";
        
        int i;
        int pos, length, row;
        String str, tmpString, chromosome;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens, tokensPrime;
        
        if (MIN_SV_LENGTH < 2*SINGLE_BREAKEND_LENGTH) {
            System.err.println("Inconsistent input: "+MIN_SV_LENGTH+"<"+(2*SINGLE_BREAKEND_LENGTH));
            System.exit(1);
        }
        currentChromosome=""; sb = new StringBuilder();
        tokensPrime = new String[10];
        br = new BufferedReader(new FileReader(INPUT_VCF));
        bw = new BufferedWriter(new FileWriter(OUTPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==VCFconstants.COMMENT) {
                if (str.substring(0,6).equals("#CHROM")) bw.write("##INFO=<ID=MATEID,Number=A,Type=String,Description=\"ID of mate breakends\">\n");
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            row=svType2Row(VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR));
            if (row==-1) {
				str=br.readLine();
				continue;
            }
            pos=Integer.parseInt(tokens[1]);
            tmpString=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (tmpString!=null) length=Integer.parseInt(tmpString);
            else length=Math.max(tokens[3].length(),tokens[4].length())-1;
            if (length<MIN_SV_LENGTH) {
				str=br.readLine();
				continue;
            }
            chromosome=tokens[0];
            pos=Integer.parseInt(tokens[1]);
            if (row==0) {  // DEL
                // Left side
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=tokens[1];
                tokensPrime[2]=tokens[2];
                tokensPrime[3]=""+getChar(chromosome,pos-1/*zero-based*/);
                tokensPrime[4]=getChar(chromosome,pos-1/*zero-based*/)+"["+chromosome+":"+(pos+length+1)+"[";
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;"+MATEID_STR+"="+NEW_MATE_PREFIX+tokens[2];
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
                // Right side
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=""+(pos+length+1);
                tokensPrime[2]=NEW_MATE_PREFIX+tokens[2];
                tokensPrime[3]=""+getChar(chromosome,pos+length+1-1/*zero-based*/);
                tokensPrime[4]="]"+chromosome+":"+pos+"]"+getChar(chromosome,pos+length+1-1/*zero-based*/);
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;"+MATEID_STR+"="+tokens[2];
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
            }
            else if (row==1) {  // INV
                // Adjacency 1, left side.
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=tokens[1];
                tokensPrime[2]=tokens[2];
                tokensPrime[3]=""+getChar(chromosome,pos-1/*zero-based*/);
                tokensPrime[4]=getChar(chromosome,pos-1/*zero-based*/)+"]"+chromosome+":"+(pos+length)+"]";
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;"+MATEID_STR+"="+NEW_MATE_PREFIX+tokens[2];
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
                // Adjacency 1, right side.
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=""+(pos+length);
                tokensPrime[2]=NEW_MATE_PREFIX+tokens[2];
                tokensPrime[3]=""+getChar(chromosome,pos+length-1/*zero-based*/);
                tokensPrime[4]=getChar(chromosome,pos+length-1/*zero-based*/)+"]"+chromosome+":"+pos+"]";
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;"+MATEID_STR+"="+tokens[2];
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
                // Adjacency 2, left side.
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=""+(pos+1);
                tokensPrime[2]=tokens[2]+SECOND_ADJACENCY_SUFFIX;
                tokensPrime[3]=""+getChar(chromosome,pos+1-1/*zero-based*/);
                tokensPrime[4]="["+chromosome+":"+(pos+length+1)+"["+getChar(chromosome,pos+1-1/*zero-based*/);
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;"+MATEID_STR+"="+NEW_MATE_PREFIX+tokens[2]+SECOND_ADJACENCY_SUFFIX;
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
                // Adjacency 1, right side.
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=""+(pos+length+1);
                tokensPrime[2]=NEW_MATE_PREFIX+tokens[2]+SECOND_ADJACENCY_SUFFIX;
                tokensPrime[3]=""+getChar(chromosome,pos+length+1-1/*zero-based*/);
                tokensPrime[4]="["+chromosome+":"+(pos+1)+"["+getChar(chromosome,pos+length+1-1/*zero-based*/);
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;"+MATEID_STR+"="+tokens[2]+SECOND_ADJACENCY_SUFFIX;
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
            }
            else if (row==2) {  // DUP
                // Left side
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=""+(pos+1);
                tokensPrime[2]=tokens[2];
                tokensPrime[3]=""+getChar(chromosome,pos+1-1/*zero-based*/);
                tokensPrime[4]="]"+chromosome+":"+(pos+length)+"]"+getChar(chromosome,pos+1-1/*zero-based*/);
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;"+MATEID_STR+"="+NEW_MATE_PREFIX+tokens[2];
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
                // Right side
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=""+(pos+length);
                tokensPrime[2]=NEW_MATE_PREFIX+tokens[2];
                tokensPrime[3]=""+getChar(chromosome,pos+length-1/*zero-based*/);
                tokensPrime[4]=getChar(chromosome,pos+length-1/*zero-based*/)+"["+chromosome+":"+(pos+1)+"[";
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;"+MATEID_STR+"="+tokens[2];
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
            }
            else if (row==3 && tokens[4].charAt(0)!='<' && tokens[4].length()-1>=MIN_SV_LENGTH) {  // INS
                // Adjacency 1
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=tokens[1];
                tokensPrime[2]=tokens[2];
                tokensPrime[3]=""+getChar(chromosome,pos-1/*zero-based*/);
                tokensPrime[4]=getChar(chromosome,pos-1/*zero-based*/)+tokens[4].substring(1,1+SINGLE_BREAKEND_LENGTH)+".";
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;";
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
                // Adjacency 2
                tokensPrime[0]=tokens[0];
                tokensPrime[1]=""+(pos+1);
                tokensPrime[2]=NEW_MATE_PREFIX+tokens[2];
                tokensPrime[3]=""+getChar(chromosome,pos+1-1/*zero-based*/);
                tokensPrime[4]="."+tokens[4].substring(tokens[4].length()-SINGLE_BREAKEND_LENGTH)+getChar(chromosome,pos+1-1/*zero-based*/);
                tokensPrime[5]=tokens[5];
                tokensPrime[6]=tokens[6];
                tokensPrime[7]="SVTYPE=BND;";
                tokensPrime[8]=tokens[8];
                tokensPrime[9]=tokens[9];
                bw.write(tokensPrime[0]);
                for (i=1; i<tokensPrime.length; i++) bw.write("\t"+tokensPrime[i]);
                bw.newLine();
            }
            str=br.readLine();
        }
        br.close(); bw.close();
    }
    
    
    /**
     * Remark: the procedure ensures that $chr$ is loaded in global variable
     * $sb$.
     */
    private static final char getChar(String chr, int first) throws IOException {
        String str;
        BufferedReader br;
        
        if (!chr.equalsIgnoreCase(currentChromosome)) {
            currentChromosome=chr;
            sb.delete(0,sb.length());
            br = new BufferedReader(new FileReader(CHROMOSOMES_DIR+"/"+chr+".fa"));
            str=br.readLine();  // Skipping FASTA header
            str=br.readLine();
            while (str!=null) {
                sb.append(str);
                str=br.readLine();
            }
            br.close();
        }
        if (first<0) first=0;
        if (first>sb.length()-1) first=sb.length()-1;
        return sb.charAt(first);
    }
    
    
	private static final int svType2Row(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(VCFconstants.DEL_STR) || 
			 type.equalsIgnoreCase(VCFconstants.DEL_ME_STR)
		   ) return 0;
		else if (type.equalsIgnoreCase(VCFconstants.INV_STR)) return 1;
        else if ( type.equalsIgnoreCase(VCFconstants.DUP_STR) ||
			      type.equalsIgnoreCase(VCFconstants.DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(VCFconstants.DUP_INT_STR) ||
                  type.equalsIgnoreCase(VCFconstants.CNV_STR)
			    ) return 2;
        else if ( type.equalsIgnoreCase(VCFconstants.INS_STR) ||
                  type.equalsIgnoreCase(VCFconstants.INS_ME_STR) ||
                  type.equalsIgnoreCase(VCFconstants.INS_NOVEL_STR)
                ) return 3;
		else return -1;
	}

}