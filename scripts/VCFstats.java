import java.io.*;


/**
 * 
 */
public class VCFstats {
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final int MAX_LENGTH = Integer.parseInt(args[1]);
        final int LENGTH_QUANTUM = Integer.parseInt(args[2]);
        
        final char COMMENT = '#';
        final String VCF_HEADER = "#CHROM";
        final String GT_SEPARATOR = ":";
        final int N_CALLERS = 3;
        final int RARE_MAX_NSAMPLES = 5;
        
        int i, j, n;
        int type, length, frequency, sampleFrequency, caller;
        int nCalls, nIns, nDel, nCalls_unsupported, nIns_unsupported, nDel_unsupported;
        String str, field;
        BufferedReader br;
        BufferedWriter bw;
        int[] supp = new int[N_CALLERS];
        int[] insHistogram, delHistogram;
        int[] insPerIndividual, delPerIndividual;  // Sites
        int[] af, sf;
        int[][] rareCalls;  // By sample freq.
        String[] tokens, tokensPrime;
        
        insHistogram = new int[1+MAX_LENGTH/LENGTH_QUANTUM];
        delHistogram = new int[1+MAX_LENGTH/LENGTH_QUANTUM];
        insPerIndividual = new int[0]; delPerIndividual = new int[0]; af = new int[0]; sf = new int[0];
        rareCalls = new int[1+RARE_MAX_NSAMPLES][N_CALLERS+1];
        nCalls=0; nIns=0; nDel=0; nCalls_unsupported=0; nIns_unsupported=0; nDel_unsupported=0;
        bw = new BufferedWriter(new FileWriter("support.txt"));
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                if (str.substring(0,VCF_HEADER.length()).equals(VCF_HEADER)) {
                    tokens=str.split("\t");
                    insPerIndividual = new int[tokens.length-9];
                    delPerIndividual = new int[tokens.length-9];
                    af = new int[1+((tokens.length-9)<<1)];
                    sf = new int[1+tokens.length-9];
                }
                str=br.readLine();
                continue;
            }
            nCalls++;
            if (nCalls%100000==0) System.err.println("Processed "+nCalls+" calls");
            tokens=str.split("\t");
            
            // SVTYPE
            type=-1;
            field=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (field!=null) type=VCFconstants.getType_infoField(field);
            else type=VCFconstants.getType_altField(tokens[4]);
            if (type==-1) {
                if (tokens[4].length()>tokens[3].length()) type=VCFconstants.TYPE_INSERTION;
                else if (tokens[3].length()>tokens[4].length()) type=VCFconstants.TYPE_DELETION;
            }
            if (type==VCFconstants.TYPE_INSERTION) nIns++;
            else if (type==VCFconstants.TYPE_DELETION) nDel++;
            
            // SVLEN
            field=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (field!=null) {
                length=Integer.parseInt(field);
                if (length>MAX_LENGTH) length=MAX_LENGTH;
                if (type==VCFconstants.TYPE_INSERTION) insHistogram[length/LENGTH_QUANTUM]++;
                else if (type==VCFconstants.TYPE_DELETION) delHistogram[length/LENGTH_QUANTUM]++;
            }
            
            // AF and SUPP
            frequency=0; sampleFrequency=0;
            for (i=0; i<N_CALLERS; i++) supp[i]=0;
            for (i=9; i<tokens.length; i++) {
                if (tokens[i].indexOf("0/1")>=0 || tokens[i].indexOf("1/0")>=0 || tokens[i].indexOf("0|1")>=0 || tokens[i].indexOf("1|0")>=0) {
                    frequency++; sampleFrequency++;
                    if (type==VCFconstants.TYPE_INSERTION) insPerIndividual[i-9]++;
                    else if (type==VCFconstants.TYPE_DELETION) delPerIndividual[i-9]++;
                }
                else if (tokens[i].indexOf("1/1")>=0 || tokens[i].indexOf("1|1")>=0) {
                    frequency+=2; sampleFrequency++;
                    if (type==VCFconstants.TYPE_INSERTION) insPerIndividual[i-9]++;
                    else if (type==VCFconstants.TYPE_DELETION) delPerIndividual[i-9]++;
                }
                tokensPrime=tokens[i].split(GT_SEPARATOR);
                n=tokensPrime.length;
                for (j=0; j<N_CALLERS; j++) {
                    if (tokensPrime[n-N_CALLERS+j].charAt(0)=='1') supp[j]++;
                }
            }
            af[frequency]++; sf[sampleFrequency]++;
            for (i=0; i<N_CALLERS; i++) bw.write(supp[i]+",");
            
            // Rare calls supported by a single caller
            if (sampleFrequency==0) {
                nCalls_unsupported++;
                if (type==VCFconstants.TYPE_INSERTION) nIns_unsupported++;
                else if (type==VCFconstants.TYPE_DELETION) nDel_unsupported++;
            }
            if (sampleFrequency<=RARE_MAX_NSAMPLES) {
                rareCalls[sampleFrequency][N_CALLERS]++;
                caller=-1;
                for (i=9; i<tokens.length; i++) {
                    tokensPrime=tokens[i].split(GT_SEPARATOR);
                    n=tokensPrime.length;
                    for (j=0; j<N_CALLERS; j++) {
                        if (tokensPrime[n-N_CALLERS+j].charAt(0)=='1') {
                            if (caller==-1) caller=j;
                            else if (caller!=j) { caller=-2; break; }
                        }
                    }
                    if (caller==-2) break;
                }
                if (caller>=0) rareCalls[sampleFrequency][caller]++;
            }
            bw.newLine();
            
            // Next iteration
            str=br.readLine();
        }
        bw.close();
        
        // Outputting
        System.out.println("nCalls="+nCalls+" nIns="+nIns+" nDel="+nDel+" other types="+(nCalls-nIns-nDel));
        System.out.println("nCalls_unsupported="+nCalls_unsupported+" nIns_unsupported="+nIns_unsupported+" nDel_unsupported="+nDel_unsupported+" otherTypes_unsupported="+(nCalls_unsupported-nIns_unsupported-nDel_unsupported));
        bw = new BufferedWriter(new FileWriter("lengthHistogram.txt"));
        for (i=0; i<insHistogram.length; i++) bw.write(insHistogram[i]+","+delHistogram[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter("individualHistogram.txt"));
        for (i=0; i<insPerIndividual.length; i++) bw.write(insPerIndividual[i]+","+delPerIndividual[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter("af.txt"));
        for (i=0; i<af.length; i++) bw.write(af[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter("sf.txt"));
        for (i=0; i<sf.length; i++) bw.write(sf[i]+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter("rareCalls.txt"));
        for (i=0; i<rareCalls.length; i++) {
            bw.write(i+"");
            for (j=0; j<rareCalls[i].length; j++) bw.write(","+rareCalls[i][j]);
            bw.newLine();
        }
        bw.close();
    }
    
}