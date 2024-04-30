import java.io.*;


/**
 * Given an inter-sample VCF with SUPP fields for every sample, the program 
 * draws the graph of all intra-sample and inter-sample merges.
 */
public class MergeGraph {
    /**
     *
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final String CHR = args[1];
        
        final int MASK_PBSV = 1;
        final int MASK_SNIFFLES = 2;
        final int MASK_PAV = 4;
        
        final String GT_SEPARATOR = ":";
        final String COLOR_PBSV = "#FF0000";
        final String COLOR_SNIFFLES = "#0000FF";
        final String COLOR_PAV = "#00FF00";
        final String COLOR_PBSV_SNIFFLES = "#FF00FF";
        final String COLOR_PBSV_PAV = "#FFFF00";
        final String COLOR_SNIFFLES_PAV = "#00FFFF";
        final String COLOR_ALL = "#000000";
        
        final double SCALE_FACTOR = 100.0;
        
        int i, j;
        int length, type, supp, suppIntersample, suppPbsv, suppSniffles, suppPav, sample, id, idGenerator;
        String str, field, color;
        BufferedReader br;
        String[] tokens, tokensPrime;
        
        System.out.println("digraph G {");
        idGenerator=-1;
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') { str=br.readLine(); continue; }
            tokens=str.split("\t");
            if (!tokens[0].equalsIgnoreCase(CHR)) break;
            id=++idGenerator;
            field=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (field!=null) length=Integer.parseInt(field);
            else length=tokens[3].length()-tokens[4].length();
            if (length<0) length=-length;
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            tokensPrime=tokens[8].split(GT_SEPARATOR);
            suppIntersample=0;
            for (i=9; i<tokens.length; i++) {
                sample=i-9;
                tokensPrime=tokens[i].split(GT_SEPARATOR);
                suppPbsv=tokensPrime[tokensPrime.length-3].charAt(0)=='1'?1:0;
                suppSniffles=tokensPrime[tokensPrime.length-2].charAt(0)=='1'?1:0;
                suppPav=tokensPrime[tokensPrime.length-1].charAt(0)=='1'?1:0;
                supp=suppPbsv|(suppSniffles<<1)|(suppPav<<2);
                if (supp==0) continue;
                suppIntersample|=supp;
                if (suppPbsv!=0) {
                    System.out.println(id+"_"+sample+"_pbsv -> "+id+"_"+sample+";");
                    System.out.println(id+"_"+sample+"_pbsv [color=\""+COLOR_PBSV+"\";Size=\""+(length/SCALE_FACTOR)+"\";type=\""+type+"\"];");
                }
                if (suppSniffles!=0) {
                    System.out.println(id+"_"+sample+"_sniffles -> "+id+"_"+sample+";");
                    System.out.println(id+"_"+sample+"_sniffles [color=\""+COLOR_SNIFFLES+"\";Size=\""+(length/SCALE_FACTOR)+"\";type=\""+type+"\"];");
                }
                if (suppPav!=0) {
                    System.out.println(id+"_"+sample+"_pav -> "+id+"_"+sample+";");
                    System.out.println(id+"_"+sample+"_pav [color=\""+COLOR_PAV+"\";Size=\""+(length/SCALE_FACTOR)+"\";type=\""+type+"\"];");
                }
                switch(supp) {
                    case MASK_PBSV: color=COLOR_PBSV; break;
                    case MASK_SNIFFLES: color=COLOR_SNIFFLES; break;
                    case MASK_PAV: color=COLOR_PAV; break;
                    case (MASK_PBSV|MASK_SNIFFLES): color=COLOR_PBSV_SNIFFLES; break;
                    case (MASK_PBSV|MASK_PAV): color=COLOR_PBSV_PAV; break;
                    case (MASK_SNIFFLES|MASK_PAV): color=COLOR_SNIFFLES_PAV; break;
                    case (MASK_PBSV|MASK_SNIFFLES|MASK_PAV): color=COLOR_ALL; break;
                    default: color="0";
                }
                System.out.println(id+"_"+sample+" [color=\""+color+"\";Size=\""+(length/SCALE_FACTOR)+"\";type=\""+type+"\"];");
                System.out.println(id+"_"+sample+" -> "+id+";");
            }
            if (suppIntersample==0) System.err.println("This call is supported by no caller in any sample?! "+str);
            else {
                switch(suppIntersample) {
                    case MASK_PBSV: color=COLOR_PBSV; break;
                    case MASK_SNIFFLES: color=COLOR_SNIFFLES; break;
                    case MASK_PAV: color=COLOR_PAV; break;
                    case (MASK_PBSV|MASK_SNIFFLES): color=COLOR_PBSV_SNIFFLES; break;
                    case (MASK_PBSV|MASK_PAV): color=COLOR_PBSV_PAV; break;
                    case (MASK_SNIFFLES|MASK_PAV): color=COLOR_SNIFFLES_PAV; break;
                    case (MASK_PBSV|MASK_SNIFFLES|MASK_PAV): color=COLOR_ALL; break;
                    default: color="0";
                }
                System.out.println(id+" [color=\""+color+"\";Size=\""+(length/SCALE_FACTOR)+"\";type=\""+type+"\"];");
            }
            str=br.readLine();
        }
        br.close();
        System.out.println("}");
    }
    
}