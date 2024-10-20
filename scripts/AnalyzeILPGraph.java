import java.util.Arrays;
import java.io.*;


public class AnalyzeILPGraph {
    /**
     * Input data
     */
    private static Triplet[] triplets;
    private static int lastTriplet;
    private static int nSamples;
    
    /**
     * Reads and read clusters. Reads wth the same ID but different samples are 
     * considered distinct.
     */
    private static int nReads;
    private static int nReadClusters;
    private static int[] readClusterID;
    private static int nSampleClusters;
    private static int[] sampleClusterID;
    private static int nHaps;
    private static int nHapClusters;
    private static int[] hapClusterID;
    

    /**
     * @param args
     * 0: all the read-haplotype edges in the CSV are assumed to be distinct;
     * 2: read-haplotype weights are taken modulo this amount before being
     * compared exactly.
     */
	public static void main(String[] args) throws IOException {
		final String INPUT_CSV = args[0];
        final String OUTPUT_DIR = args[1];
        final int WEIGHT_QUANTUM = Integer.parseInt(args[2]);
		
        final String BASENAME = INPUT_CSV.substring(INPUT_CSV.lastIndexOf("/")+1,INPUT_CSV.lastIndexOf("."));
        final int CAPACITY = 100000;  // Arbitrary
        
		String str;
		BufferedReader br;
        String[] tokens;
        
        System.out.println("Loading CSV file...");
        triplets = new Triplet[CAPACITY];
        br = new BufferedReader(new FileReader(INPUT_CSV));
        lastTriplet=-1;
		str=br.readLine(); str=br.readLine();  // Skipping header
        while (str!=null) {
            tokens=str.split(",");
            lastTriplet++;
            if (lastTriplet==triplets.length) {
                Triplet[] newArray = new Triplet[triplets.length+CAPACITY];
                System.arraycopy(triplets,0,newArray,0,triplets.length);
                triplets=newArray;
            }
            triplets[lastTriplet] = new Triplet(tokens[0],tokens[1],tokens[3],Integer.parseInt(tokens[5])/WEIGHT_QUANTUM);
            str=br.readLine();
        }
	    br.close();

        // Computing clusters of identical reads, haplotypes, samples.
        System.out.println("Finding identical reads...");
        identicalReads();
        System.out.println("Finding identical samples...");
        identicalSamples();
        System.out.println("Finding identical haplotypes...");
        identicalHaplotypes();
        System.err.println(nReads+"\t"+nReadClusters+"\t"+(((double)nReads)/nReadClusters)+"\t"+nHaps+"\t"+nHapClusters+"\t"+(((double)nHaps)/nHapClusters)+"\t"+nSamples+"\t"+nSampleClusters+"\t"+(((double)nSamples)/nSampleClusters));
	}
    
    
    /**
     * Two reads are considered identical iff they connect to the same
     * haplotypes with the same (possibly quantized) weights.
     *
     * Remark: the procedure sorts global array $triplets$ by (sample,read,
     * haplotype), and sets global variables $nSamples,clusterID,nClusters$.
     * $clusterID$ is set to the cluster ID of each read in the order in which 
     * it appears in the sorted $triplets$.
     */
    private static final void identicalReads() throws IOException {
        final int CAPACITY = 20;  // Even
        
        boolean found;
        int i, j, k;
        String lastSample, lastRead;
        int[] lastNeighbor;
        String[][] neighbors;
        
        Triplet.order=Triplet.ORDER_SAMPLE_READ_HAP;
        Arrays.sort(triplets,0,lastTriplet+1);
        j=-1; nReads=0;
        for (i=0; i<=lastTriplet; i++) {
            if (j==-1 || !triplets[i].sample.equals(triplets[j].sample) || !triplets[i].read.equals(triplets[j].read)) { nReads++; j=i; }
        }
        
        // Collecting all read-haplotype edges
        neighbors = new String[nReads][CAPACITY];
        lastNeighbor = new int[nReads];
        for (i=0; i<lastNeighbor.length; i++) lastNeighbor[i]=-1;
        j=0; 
        lastSample=triplets[0].sample; nSamples=1; 
        lastRead=triplets[0].read;
        neighbors[0][0]=triplets[0].hap; neighbors[0][1]=triplets[0].count+"";
        lastNeighbor[0]=1;
        for (i=1; i<=lastTriplet; i++) {
            if (triplets[i].sample.equals(lastSample) && triplets[i].read.equals(lastRead)) {
                if (lastNeighbor[j]+2>=neighbors[j].length) {
                    String[] newArray = new String[neighbors[j].length+CAPACITY];
                    System.arraycopy(neighbors[j],0,newArray,0,neighbors[j].length);
                    neighbors[j]=newArray;
                }
                neighbors[j][++lastNeighbor[j]]=triplets[i].hap;
                neighbors[j][++lastNeighbor[j]]=triplets[i].count+"";
            }
            else {
                if (!triplets[i].sample.equals(lastSample)) { nSamples++; lastSample=triplets[i].sample; }
                lastRead=triplets[i].read;
                j++;
                lastNeighbor[j]=-1;
                neighbors[j][++lastNeighbor[j]]=triplets[i].hap;
                neighbors[j][++lastNeighbor[j]]=triplets[i].count+"";
            }
        }
        
        // Clustering reads
        readClusterID = new int[nReads];
        for (i=0; i<readClusterID.length; i++) readClusterID[i]=-1;
        nReadClusters=0;
        for (i=0; i<nReads; i++) {
            if (readClusterID[i]!=-1) continue;
            readClusterID[i]=nReadClusters++;
            for (j=i+1; j<nReads; j++) {
                if (readClusterID[j]!=-1 || lastNeighbor[j]!=lastNeighbor[i]) continue;
                found=true;
                for (k=0; k<=lastNeighbor[i]; k+=2) {
                    if (!neighbors[i][k].equals(neighbors[j][k]) || Integer.parseInt(neighbors[i][k+1])!=Integer.parseInt(neighbors[j][k+1])) { found=false; break; }
                }
                if (found) readClusterID[j]=readClusterID[i];
            }
        }
    }
    
    
    /**
     * Two samples are considered identical iff their reads belong to the same
     * set of read clusters.
     *
     * Remark: the procedure assumes that global array $triplets$ is sorted by
     * (sample,read,hap), and that global variable $readClusterID$ has already
     * been initialized.
     *
     * Remark: the procedure sets global variables $sampleClusterID,
     * nSampleClusters$.
     */
    private static final void identicalSamples() {
        boolean found;
        int i, j, x;
        String lastSample, lastRead;
        boolean[][] matrix;
        
        // Building the sample-cluster matrix
        matrix = new boolean[nSamples][nReadClusters];
        lastSample=""; lastRead=""; i=-1; j=-1;
        for (x=0; x<lastTriplet; x++) {
            if (!triplets[x].sample.equals(lastSample)) { 
                i++; lastSample=triplets[x].sample; 
                j++; lastRead=triplets[x].read;
            }
            else if (!triplets[x].read.equals(lastRead)) { j++; lastRead=triplets[x].read; } 
            matrix[i][readClusterID[j]]=true;
        }
        
        // Counting identical samples
        sampleClusterID = new int[nSamples];
        for (i=0; i<nSamples; i++) sampleClusterID[i]=-1;
        nSampleClusters=0;
        for (i=0; i<nSamples; i++) {
            if (sampleClusterID[i]!=-1) continue;
            sampleClusterID[i]=nSampleClusters++;
            for (j=i+1; j<nSamples; j++) {
                found=true;
                for (x=0; x<nReadClusters; x++) {
                    if (matrix[i][x]!=matrix[j][x]) { found=false; break; }
                }
                if (found) sampleClusterID[j]=sampleClusterID[i];
            }
        }
    }
    
    
    /**
     * Two haplotypes are considered identical iff they connect to the same
     * reads with the same (possibly quantized) weights.
     *
     * Remark: the procedure sorts global variable $triplets$ by (hap,sample,
     * read) and sets global variables $nHaps,hapClusterID,nHapClusters$.
     * 
     * Remark: every (read,hap) pair is assumed to occur just once in 
     * $triplets$.
     */
    private static final void identicalHaplotypes() throws IOException {
        final int CAPACITY = 30;  // Multiple of 3
        
        boolean found;
        int i, j, k;
        int nClusters;
        String lastHaplotype;
        int[] lastNeighbor;
        String[][] neighbors;
        
        Triplet.order=Triplet.ORDER_HAP_SAMPLE_READ;
        Arrays.sort(triplets,0,lastTriplet+1);
        j=-1; nHaps=0;
        for (i=0; i<=lastTriplet; i++) {
            if (j==-1 || !triplets[i].hap.equals(triplets[j].hap)) { nHaps++; j=i; }
        }
        
        // Collecting all haplotype-read edges
        neighbors = new String[nHaps][CAPACITY];
        lastNeighbor = new int[nHaps];
        for (i=0; i<lastNeighbor.length; i++) lastNeighbor[i]=-1;
        j=0; lastHaplotype=triplets[0].hap;
        neighbors[0][0]=triplets[0].sample; neighbors[0][1]=triplets[0].read; neighbors[0][2]=triplets[0].count+"";
        lastNeighbor[0]=2;
        for (i=1; i<=lastTriplet; i++) {
            if (triplets[i].hap.equals(lastHaplotype)) {
                if (lastNeighbor[j]+3>=neighbors[j].length) {
                    String[] newArray = new String[neighbors[j].length+CAPACITY];
                    System.arraycopy(neighbors[j],0,newArray,0,neighbors[j].length);
                    neighbors[j]=newArray;
                }
                neighbors[j][++lastNeighbor[j]]=triplets[i].sample;
                neighbors[j][++lastNeighbor[j]]=triplets[i].read;
                neighbors[j][++lastNeighbor[j]]=triplets[i].count+"";
            }
            else {
                lastHaplotype=triplets[i].hap;
                j++;
                lastNeighbor[j]=-1;
                neighbors[j][++lastNeighbor[j]]=triplets[i].sample;
                neighbors[j][++lastNeighbor[j]]=triplets[i].read;
                neighbors[j][++lastNeighbor[j]]=triplets[i].count+"";
            }
        }
        
        // Clustering haplotypes
        hapClusterID = new int[nHaps];
        for (i=0; i<hapClusterID.length; i++) hapClusterID[i]=-1;
        nHapClusters=0;
        for (i=0; i<nHaps; i++) {
            if (hapClusterID[i]!=-1) continue;
            hapClusterID[i]=nHapClusters++;
            for (j=i+1; j<nHaps; j++) {
                if (hapClusterID[j]!=-1 || lastNeighbor[j]!=lastNeighbor[i]) continue;
                found=true;
                for (k=0; k<=lastNeighbor[i]; k+=3) {
                    if (!neighbors[i][k].equals(neighbors[j][k]) || !neighbors[i][k+1].equals(neighbors[j][k+1]) || Integer.parseInt(neighbors[i][k+2])!=Integer.parseInt(neighbors[j][k+2])) { found=false; break; }
                }
                if (found) hapClusterID[j]=hapClusterID[i];
            }
        }
    }
    
    
    public static class Triplet implements Comparable {
        public final static int ORDER_SAMPLE_READ_HAP = 0;
        public final static int ORDER_HAP_SAMPLE_READ = 1;
        public static int order;
        
        public String sample, read, hap;
        public int count;
        
        public Triplet(String s, String r, String h, int c) {
            this.sample=s; this.read=r; this.hap=h; this.count=c;
        }
        
        public int compareTo(Object other) {
            Triplet otherTriplet = (Triplet)other;
            
            final int s = sample.compareTo(otherTriplet.sample);    
            final int r = read.compareTo(otherTriplet.read);
            final int h = hap.compareTo(otherTriplet.hap);
            
            if (order==ORDER_SAMPLE_READ_HAP) {
                if (s<0) return -1;
                else if (s>0) return 1;
                if (r<0) return -1;
                else if (r>0) return 1;
                if (h<0) return -1;
                else if (h>0) return 1;
                else return 0;
            }
            else if (order==ORDER_HAP_SAMPLE_READ) {
                if (h<0) return -1;
                else if (h>0) return 1;
                if (s<0) return -1;
                else if (s>0) return 1;
                if (r<0) return -1;
                else if (r>0) return 1;
                else return 0;
            }
            else return 0;
        }
        
        public boolean equals(Object other) {
            Triplet otherTriplet = (Triplet)other;
            return sample.equalsIgnoreCase(otherTriplet.sample) && read.equalsIgnoreCase(otherTriplet.read) && hap.equalsIgnoreCase(otherTriplet.hap) && count==otherTriplet.count;
        }
        
        public String toString() {
            return sample+","+read+","+hap+" ("+count+")";
        }
    }
    
    
    
    
    
    
    
    /*
    private static final void printGraph(Triplet[] readHaplotype, int last_readHaplotype, int last_haplotypeSample, String outputDir, String basename) throws IOException {
        int i;
        BufferedWriter bwAll, bwReads, bwSamples;
        
        
            tokens[1]=tokens[1].replace('/','_');  // For DOT
            tokens[3]=tokens[3].replace('>','_');  // For DOT
            tokens[3]=tokens[3].replace('<','-');  // For DOT
        
        
        
            System.out.println("Printing the graph...");
            printReadHaplotype(readHaplotype,last_readHaplotype);
            
            // Computing distinct haplotype-sample edges
            Triplet.order=Triplet.ORDER_HAP_SAMPLE_READ;
            Arrays.sort(readHaplotype,0,last_readHaplotype+1);
            j=0;
            for (i=1; i<=last_readHaplotype; i++) {
                if (readHaplotype[i].hap.equals(readHaplotype[j].hap) && readHaplotype[i].sample.equals(readHaplotype[j].sample)) readHaplotype[j].count+=readHaplotype[i].count;
                else {
                    j++;
                    readHaplotype[j].sample=readHaplotype[i].sample;
                    readHaplotype[j].read=readHaplotype[i].read;
                    readHaplotype[j].hap=readHaplotype[i].hap;
                }
            }
            last_readHaplotype=j;
            
            printHaplotypeSample(readHaplotype,last_readHaplotype);
        
        
        
        
        
        
        
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
    */
    

}