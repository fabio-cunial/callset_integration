import java.io.*;


/**
 * Transforms a VCF into the TSV needed by sv-merger (printed to STDOUT):
 *
 * chr \t start \t end \t svID \t sampleID \t caller \t svType \t svLength
 *
 */
public class VCF2SVMerger {

    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final String SAMPLE_ID = args[1];
        final String CALLER_ID = args[2];
        
        int type, pos, end, length, nCalls, nDiscarded_chromosome, nDiscarded_type, nDiscarded_length;
        String str, field;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine(); nCalls=0; nDiscarded_chromosome=0; nDiscarded_type=0; nDiscarded_length=0;
        while (str!=null) {
            if (str.charAt(0)==VCFconstants.COMMENT) {
                str=br.readLine();
                continue;
            }
            if (nCalls%10000==0) System.err.println("Processed "+nCalls+" calls");
            nCalls++;
            tokens=str.split("\t");
            if (!VCFconstants.isChromosome(tokens[0])) {
                nDiscarded_chromosome++;
                str=br.readLine();
                continue;
            }
            
            // Begin site
            pos=Integer.parseInt(tokens[1]);
            
            // SV type
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type!=VCFconstants.TYPE_INSERTION && type!=VCFconstants.TYPE_DELETION && type!=VCFconstants.TYPE_DUPLICATION && type!=VCFconstants.TYPE_INVERSION) {
                nDiscarded_type++;
                str=br.readLine();
                continue;
            }
            
            // End site and SV length
            field=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (field!=null) {
                length=Integer.parseInt(field);
                if (length<0) length=-length;
                end=pos+length;  // This is necessary also for INS
            }
            else {
                if (tokens[4].charAt(0)!='<' && tokens[4].charAt(tokens[4].length()-1)!='>') {
                    length=tokens[4].length()-tokens[3].length();
                    if (length<0) length=-length;
                    end=pos+length;
                }
                else {
                    field=VCFconstants.getField(tokens[7],VCFconstants.END_STR);
                    if (field!=null) {
                        end=Integer.parseInt(field);
                        length=end-pos;
                    }
                    else {
                        nDiscarded_length++;
                        str=br.readLine();
                        continue;
                    }
                }
            }
            
            System.out.println(tokens[0]+"\t"+pos+"\t"+end+"\t"+CALLER_ID+"-"+tokens[2]+"\t"+SAMPLE_ID+"\t"+CALLER_ID+"\t"+VCFconstants.type2str(type)+"\t"+length);
            str=br.readLine();
        }
        br.close();
        System.err.println("Discarded "+nDiscarded_chromosome+" calls on a non-standard chromosome ("+((100.0*nDiscarded_chromosome)/nCalls)+"%)");
        System.err.println("Discarded "+nDiscarded_type+" calls of type different from INS,DEL ("+((100.0*nDiscarded_type)/nCalls)+"%)");
        System.err.println("Discarded "+nDiscarded_length+" calls with unknown length/end ("+((100.0*nDiscarded_length)/nCalls)+"%).");
    }
    
}