import java.io.*;


/**
 * Given a VCF that contains just BNDs, the program:
 * - Removes virtual telomeric breakends.
 * - Ensures that REF and ALT use the same first character, and that such a
 *   character is the same as in the ref.
 * - Ensures that every non-single BND has a mate.
 */
public class CleanBNDs {
    /**
     * Chromosome buffer
     */
    private static String CHROMOSOMES_DIR;
    private static String currentChromosome;  // ID
    private static StringBuilder sb;  // Sequence
    
    private static String[] chromosomeIDs;
    private static int[] chromosomeLengths;
    private static int nChromosomes;
    
    
    /**
     * @param args 
     * 1=directory that contains one $Z.fa$ file for every chromosome $Z$. These
     * files can be created from a single reference file by doing e.g.: 
     *
     * cat ref.fasta | awk '{ if (substr($0,1,1)==">") { filename=(substr($1,2) ".fa") } print $0 >> filename; close(filename) }'
     *
     * The directory is also assumed to contain a file $index.fai$.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        CHROMOSOMES_DIR=args[1];
        final String OUTPUT_VCF = args[2];
        
        final char COMMENT = '#';
        final String MATEID_STR = "MATEID";
        final String NEW_MATE_PREFIX = "mate_of_";
        boolean isSingle, hasMate;
        int pos, pos2, length, nCalls;
        String chromosome, chromosome2, id, ref, alt, qual, filter, info, format, gt;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        
        loadChromosomeLengths(CHROMOSOMES_DIR+"/index.fai");
        
        // Fixing the VCF
        nCalls=0; currentChromosome=""; sb = new StringBuilder();
        br = new BufferedReader(new FileReader(INPUT_VCF));
        bw = new BufferedWriter(new FileWriter(OUTPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                if (str.substring(0,6).equals("#CHROM")) bw.write("##INFO=<ID=MATEID,Number=A,Type=String,Description=\"ID of mate breakends\">\n");
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            nCalls++;
            tokens=str.split("\t");
            chromosome=tokens[0];
            pos=Integer.parseInt(tokens[1]);
            if (pos==0 || pos==getChromosomeLength(chromosome)+1) {
                // Virtual telomeric breakend
                str=br.readLine();
                continue;
            }
            id=tokens[2]; ref=tokens[3]; alt=tokens[4]; qual=tokens[5]; filter=tokens[6]; info=tokens[7]; format=tokens[8]; gt=tokens[9];
            isSingle=alt.charAt(0)=='.'||alt.charAt(alt.length()-1)=='.';
            hasMate=info.indexOf(MATEID_STR)>=0;
            if (!isSingle) {
                chromosome2=alt2chromosome(alt);
                pos2=alt2pos(alt);
                if (pos2==0 || pos2==getChromosomeLength(chromosome2)+1) {
                    // Virtual telomeric breakend
                    str=br.readLine();
                    continue;
                }
            }
            else { chromosome2=null; pos2=0; }
            
            // Fixing REF and ALT, if needed.
            if (ref.charAt(0)=='N' || ref.charAt(0)=='n') ref=getChar(chromosome,pos-1/*zero-based*/)+ref.substring(1);
            length=alt.length();
            if (alt.charAt(0)=='N' || alt.charAt(0)=='n') alt=getChar(chromosome,pos-1/*zero-based*/)+alt.substring(1);
            if (alt.charAt(length-1)=='N' || alt.charAt(length-1)=='n') alt=alt.substring(0,length-1)+getChar(chromosome,pos-1/*zero-based*/);
            if (!isSingle && !hasMate) info+=";"+MATEID_STR+"="+NEW_MATE_PREFIX+id;
            bw.write(chromosome); bw.write("\t"+pos); bw.write("\t"+id); bw.write("\t"+ref); bw.write("\t"+alt); bw.write("\t"+qual); bw.write("\t"+filter); bw.write("\t"+info); bw.write("\t"+format); bw.write("\t"+gt); bw.newLine();
            
            // Writing the new mate record, if any.
            if (!isSingle && !hasMate) { bw.write(chromosome2); bw.write("\t"+pos2); bw.write("\t"+NEW_MATE_PREFIX+id); bw.write("\t"+getChar(chromosome2,pos2-1/*zero-based*/)); bw.write("\t"+getMateAlt(alt,chromosome,pos,getChar(chromosome2,pos2-1/*zero-based*/))); bw.write("\t"+qual); bw.write("\t"+filter); bw.write("\tSVTYPE=BND;"+MATEID_STR+"="+id); bw.write("\t"+format); bw.write("\t"+gt); bw.newLine(); }
            
            // Next record
            str=br.readLine();
        }
        br.close(); bw.close();
    }
    
    
    private static final void loadChromosomeLengths(String input) throws IOException {
        int i;
        String str;
        BufferedReader br;
        String[] tokens;
        
        nChromosomes=0;
        br = new BufferedReader(new FileReader(input));
        str=br.readLine();
        while (str!=null) { nChromosomes++; str=br.readLine(); }
        br.close();
        chromosomeIDs = new String[nChromosomes];
        chromosomeLengths = new int[nChromosomes];
        br = new BufferedReader(new FileReader(input));
        str=br.readLine(); i=0;
        while (str!=null) {
            tokens=str.split("\t");
            chromosomeIDs[i]=tokens[0];
            chromosomeLengths[i]=Integer.parseInt(tokens[1]);
            i++;
            str=br.readLine();
        }
        br.close();
    }
    
    
    private static final int getChromosomeLength(String chr) {
        for (int i=0; i<nChromosomes; i++) {
            if (chromosomeIDs[i].equalsIgnoreCase(chr)) return chromosomeLengths[i];
        }
        return 0;
    }
    
    
    private static final String alt2chromosome(String alt) {
        int p;
        
        p=alt.indexOf("[");
        if (p<0) p=alt.indexOf("]");
        p++;
        return alt.substring(p,alt.indexOf(":",p+1));
    }


    private static final int alt2pos(String alt) {
        int p, q;
        
        p=alt.indexOf(":")+1;
        q=alt.indexOf("[",p+1);
        if (q<0) q=alt.indexOf("]",p+1);
        return Integer.parseInt(alt.substring(p,q));
    }
    
    
    private static final String getMateAlt(String alt, String chromosome1, int pos1, char char2) {
        final char c = alt.charAt(0);
        
        if (c=='[') return "["+chromosome1+":"+pos1+"["+char2;
        else if (c==']') return char2+"["+chromosome1+":"+pos1+"[";
        else if (alt.indexOf("[")>=0) return "]"+chromosome1+":"+pos1+"]"+char2;
        else return char2+"]"+chromosome1+":"+pos1+"]";
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

}