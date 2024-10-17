import java.util.Arrays;
import java.io.*;


public class DrawILPGraph {

    /**
     * @param args
     * 0: all the read-haplotype edges in the CSV are assumed to be distinct;
     * 2: read-haplotype weights are divided by this amount before being
     * compared exactly.
     */
	public static void main(String[] args) throws IOException {
		final String INPUT_CSV = args[0];
        final String OUTPUT_DIR = args[1];
        final int WEIGHT_QUANTUM = Integer.parseInt(args[2]);
		
        final String BASENAME = INPUT_CSV.substring(INPUT_CSV.lastIndexOf("/")+1,INPUT_CSV.lastIndexOf("."));
        final int CAPACITY = 100000;  // Arbitrary
        
        int i, j;
        int last_readHaplotype, last_haplotypeSample;
		String str;
		BufferedReader br;
        String[] tokens;
        Pair[] readHaplotype, haplotypeSample;
        
        System.out.println("Loading the file...");
        readHaplotype = new Pair[CAPACITY];
        haplotypeSample = new Pair[CAPACITY];
        br = new BufferedReader(new FileReader(INPUT_CSV));
        last_readHaplotype=-1; last_haplotypeSample=-1;
		str=br.readLine(); str=br.readLine();  // Skipping header
        while (str!=null) {
            tokens=str.split(",");
            tokens[1]=tokens[1].replace('/','_');  // For DOT
            tokens[3]=tokens[3].replace('>','_');  // For DOT
            tokens[3]=tokens[3].replace('<','-');  // For DOT
            last_readHaplotype++;
            if (last_readHaplotype==readHaplotype.length) {
                Pair[] newArray = new Pair[readHaplotype.length+CAPACITY];
                System.arraycopy(readHaplotype,0,newArray,0,readHaplotype.length);
                readHaplotype=newArray;
            }
            readHaplotype[last_readHaplotype] = new Pair(tokens[1],tokens[3],Integer.parseInt(tokens[5])/WEIGHT_QUANTUM);
            last_haplotypeSample++;
            if (last_haplotypeSample==haplotypeSample.length) {
                Pair[] newArray = new Pair[haplotypeSample.length+CAPACITY];
                System.arraycopy(haplotypeSample,0,newArray,0,haplotypeSample.length);
                haplotypeSample=newArray;
            }
            haplotypeSample[last_haplotypeSample] = new Pair(tokens[3],tokens[0],1);          
            str=br.readLine();
        }
	    br.close();
        
        // Computing distinct haplotype-sample edges
        Arrays.sort(haplotypeSample,0,last_haplotypeSample+1);
        j=0;
        for (i=1; i<=last_haplotypeSample; i++) {
            if (haplotypeSample[i].equals(haplotypeSample[j])) haplotypeSample[j].count+=haplotypeSample[i].count;
            else {
                j++;
                haplotypeSample[j].from=haplotypeSample[i].from;
                haplotypeSample[j].to=haplotypeSample[i].to;
                haplotypeSample[j].count=haplotypeSample[i].count;
            }
        }
        last_haplotypeSample=j;
        
//System.out.println("Printing the graph...");        //printGraph(readHaplotype,last_readHaplotype,haplotypeSample,last_haplotypeSample,OUTPUT_DIR,BASENAME);
        
        // Computing clusters of identical reads and haplotypes
        identicalReads(readHaplotype,last_readHaplotype);
        System.err.print(",");
        identicalHaplotypes(readHaplotype,last_readHaplotype);
        System.err.println();
//Pair.order=Pair.ORDER_TO_FROM;
//Arrays.sort(readHaplotype,0,last_readHaplotype+1);        //identicalHaplotypes(readHaplotype,last_readHaplotype,WEIGHT_THRESHOLD);
	}
    
    
    /**
     * Two reads are considered identical iff they connect to the same
     * haplotypes with the same weights.
     */
    private static final void identicalReads(Pair[] readHaplotype, int last_readHaplotype) throws IOException {
        final int CAPACITY = 20;  // Even
        
        boolean found;
        int i, j, k;
        int nClusters, nReads;
        String lastRead;
        int[] lastNeighbor, count;
        String[][] neighbors;
        
        Pair.order=Pair.ORDER_FROM_TO;
        Arrays.sort(readHaplotype,0,last_readHaplotype+1);
        j=-1; nReads=0;
        for (i=0; i<=last_readHaplotype; i++) {
            if (j==-1 || !readHaplotype[i].from.equals(readHaplotype[j].from)) { nReads++; j=i; }
        }
        
        // Collecting all read-haplotype edges
        neighbors = new String[nReads][CAPACITY];
        lastNeighbor = new int[nReads];
        for (i=0; i<lastNeighbor.length; i++) lastNeighbor[i]=-1;
        j=0; lastRead=readHaplotype[0].from;
        neighbors[0][0]=readHaplotype[0].to; neighbors[0][1]=readHaplotype[0].count+"";
        lastNeighbor[0]=1;
        for (i=1; i<=last_readHaplotype; i++) {
            if (readHaplotype[i].from.equals(lastRead)) {
                if (lastNeighbor[j]+2>=neighbors[j].length) {
                    String[] newArray = new String[neighbors[j].length+CAPACITY];
                    System.arraycopy(neighbors[j],0,newArray,0,neighbors[j].length);
                    neighbors[j]=newArray;
                }
                neighbors[j][++lastNeighbor[j]]=readHaplotype[i].to;
                neighbors[j][++lastNeighbor[j]]=readHaplotype[i].count+"";
            }
            else {
                lastRead=readHaplotype[i].from;
                j++;
                lastNeighbor[j]=-1;
                neighbors[j][++lastNeighbor[j]]=readHaplotype[i].to;
                neighbors[j][++lastNeighbor[j]]=readHaplotype[i].count+"";
            }
        }
        
        // Clustering reads
        count = new int[nReads];
        for (i=0; i<count.length; i++) count[i]=0;
        nClusters=0;
        for (i=0; i<nReads; i++) {
            if (count[i]!=0) continue;
            count[i]=1; nClusters++;
            for (j=i+1; j<nReads; j++) {
                if (count[j]!=0 || lastNeighbor[j]!=lastNeighbor[i]) continue;
                found=true;
                for (k=0; k<=lastNeighbor[i]; k+=2) {
                    if (!neighbors[i][k].equals(neighbors[j][k]) || Integer.parseInt(neighbors[i][k+1])!=Integer.parseInt(neighbors[j][k+1])) { found=false; break; }
                }
                if (found) { count[i]++; count[j]=-1; }
            }
        }
        System.err.print(nReads+","+nClusters+","+(((double)nReads)/nClusters));
    }
    
    
    /**
     * Two haplotypes are considered identical iff they connect to the same
     * reads with the same weights.
     */
    private static final void identicalHaplotypes(Pair[] readHaplotype, int last_readHaplotype) throws IOException {
        final int CAPACITY = 20;  // Even
        
        boolean found;
        int i, j, k;
        int nClusters, nHaplotypes;
        String lastHaplotype;
        int[] lastNeighbor, count;
        String[][] neighbors;
        
        Pair.order=Pair.ORDER_TO_FROM;
        Arrays.sort(readHaplotype,0,last_readHaplotype+1);
        j=-1; nHaplotypes=0;
        for (i=0; i<=last_readHaplotype; i++) {
            if (j==-1 || !readHaplotype[i].to.equals(readHaplotype[j].to)) { nHaplotypes++; j=i; }
        }
        
        // Collecting all haplotype-read edges
        neighbors = new String[nHaplotypes][CAPACITY];
        lastNeighbor = new int[nHaplotypes];
        for (i=0; i<lastNeighbor.length; i++) lastNeighbor[i]=-1;
        j=0; lastHaplotype=readHaplotype[0].to;
        neighbors[0][0]=readHaplotype[0].from; neighbors[0][1]=readHaplotype[0].count+"";
        lastNeighbor[0]=1;
        for (i=1; i<=last_readHaplotype; i++) {
            if (readHaplotype[i].to.equals(lastHaplotype)) {
                if (lastNeighbor[j]+2>=neighbors[j].length) {
                    String[] newArray = new String[neighbors[j].length+CAPACITY];
                    System.arraycopy(neighbors[j],0,newArray,0,neighbors[j].length);
                    neighbors[j]=newArray;
                }
                neighbors[j][++lastNeighbor[j]]=readHaplotype[i].from;
                neighbors[j][++lastNeighbor[j]]=readHaplotype[i].count+"";
            }
            else {
                lastHaplotype=readHaplotype[i].to;
                j++;
                lastNeighbor[j]=-1;
                neighbors[j][++lastNeighbor[j]]=readHaplotype[i].from;
                neighbors[j][++lastNeighbor[j]]=readHaplotype[i].count+"";
            }
        }
        
        // Clustering haplotypes
        count = new int[nHaplotypes];
        for (i=0; i<count.length; i++) count[i]=0;
        nClusters=0;
        for (i=0; i<nHaplotypes; i++) {
            if (count[i]!=0) continue;
            count[i]=1; nClusters++;
            for (j=i+1; j<nHaplotypes; j++) {
                if (count[j]!=0 || lastNeighbor[j]!=lastNeighbor[i]) continue;
                found=true;
                for (k=0; k<=lastNeighbor[i]; k+=2) {
                    if (!neighbors[i][k].equals(neighbors[j][k]) || Integer.parseInt(neighbors[i][k+1])!=Integer.parseInt(neighbors[j][k+1])) { found=false; break; }
                }
                if (found) { count[i]++; count[j]=-1; }
            }
        }
        System.err.print(nHaplotypes+","+nClusters+","+(((double)nHaplotypes)/nClusters));
    }
    
    
    private static final void printGraph(Pair[] readHaplotype, int last_readHaplotype, Pair[] haplotypeSample, int last_haplotypeSample, String outputDir, String basename) throws IOException {
        int i;
        BufferedWriter bwAll, bwReads, bwSamples;
        
        bwAll = new BufferedWriter(new FileWriter(outputDir+"/"+basename+"_all.dot"));
        bwReads = new BufferedWriter(new FileWriter(outputDir+"/"+basename+"_reads.dot"));
        bwSamples = new BufferedWriter(new FileWriter(outputDir+"/"+basename+"_samples.dot"));
        bwAll.write("graph G {\n"); bwReads.write("graph G {\n"); bwSamples.write("graph G {\n");
        //bwAll.write("overlap = false;"); bwReads.write("overlap = false;"); bwSamples.write("overlap = false;");
        // Reads
        for (i=0; i<=last_readHaplotype; i++) {
            bwAll.write(readHaplotype[i].from+" [shape=point,style=filled,color=red,layer=\"reads\"];\n");
            bwReads.write(readHaplotype[i].from+" [shape=point,style=filled,color=red,layer=\"reads\"];\n");
        }
        // Haplotypes
        for (i=0; i<=last_haplotypeSample; i++) {
            bwAll.write(haplotypeSample[i].from+" [shape=point,style=filled,color=green,width=1,layer=\"haplotypes\"];\n");
            bwSamples.write(haplotypeSample[i].from+" [shape=point,style=filled,color=green,width=1,layer=\"haplotypes\"];\n");
            bwReads.write(haplotypeSample[i].to+" [shape=point,style=filled,color=green,width=1,layer=\"haplotypes\"];\n");
        }
        // Samples
        for (i=0; i<=last_haplotypeSample; i++) {
            bwAll.write(haplotypeSample[i].to+" [shape=point,style=filled,color=blue,layer=\"samples\"];\n");
            bwSamples.write(haplotypeSample[i].to+" [shape=point,style=filled,color=blue,layer=\"samples\"];\n");
        }
        // Edges
        for (i=0; i<=last_readHaplotype; i++) {
            bwAll.write(readHaplotype[i].from+" -- "+readHaplotype[i].to+" [weight="+readHaplotype[i].count+",color=red];\n");
            bwReads.write(readHaplotype[i].from+" -- "+readHaplotype[i].to+" [weight="+readHaplotype[i].count+",color=red];\n");
        }
        for (i=0; i<=last_haplotypeSample; i++) {
            bwAll.write(haplotypeSample[i].from+" -- "+haplotypeSample[i].to+" [weight="+haplotypeSample[i].count+",color=blue];\n");
            bwSamples.write(haplotypeSample[i].from+" -- "+haplotypeSample[i].to+" [weight="+haplotypeSample[i].count+",color=blue];\n");
        }
        bwAll.write("}"); bwReads.write("}"); bwSamples.write("}");
        bwAll.close(); bwReads.close(); bwSamples.close();
    }
    
    
    public static class Pair implements Comparable {
        public final static int ORDER_FROM_TO = 0;
        public final static int ORDER_TO_FROM = 1;
        
        public static int order;
        public String from, to;
        public int count;
        
        public Pair(String f, String t,int c) {
            this.from=f; this.to=t;
            this.count=c;
        }
        
        public int compareTo(Object other) {
            Pair otherPair = (Pair)other;

            final int f = from.compareTo(otherPair.from);
            final int t = to.compareTo(otherPair.to);
            if (order==ORDER_FROM_TO) {
                if (f<0) return -1;
                else if (f==0) {
                    if (t<0) return -1;
                    else if (t==0) return 0;
                    else return 1;
                }
                else return 1;
            }
            else if (order==ORDER_TO_FROM) {
                if (t<0) return -1;
                else if (t==0) {
                    if (f<0) return -1;
                    else if (f==0) return 0;
                    else return 1;
                }
                else return 1;
            }
            else return 0;
        }
        
        public boolean equals(Object other) {
            Pair otherPair = (Pair)other;
            return from.equalsIgnoreCase(otherPair.from) && to.equalsIgnoreCase(otherPair.to);
        }
        
        public String toString() {
            return from+" -> "+to+" ("+count+")";
        }
    }

}