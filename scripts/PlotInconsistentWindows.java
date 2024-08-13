import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class PlotInconsistentWindows {
    private static final int DEL = 0;
    private static final int INV = 1;
    private static final int DUP = 2;
    private static final int INS = 3;
    private static final int SNP = 4;
    private static final int REPLACEMENT = 5;
    
    /**
     * Temporary variables to store the graph of a window
     */
    private static int nNodes, nEdges;
    private static int[] lastNeighbor, color, queue;
    private static int[][] neighbors;
    
    /**
     * Temporary variables for computing independent sets
     */
    private static Interval[] sorted;
    private static int[] node2sorted;
    
    
    /**
     * @param args
     * 0: input VCF, assumed to be sorted; '.' in the GT is assumed to be equal 
     *    to '0';
     * 1: one-based sample column ID (1=first sample);
     * 2: for each VCF record: 
     *    isSV, chrID, first, last, gtId
     * 3: for each window: 
     *    firstRecord (0-based), lastRecord (0-based, incl.), nCollisions
     * 4: a graphical representation of the overlap graph of every window with 
     *    more than one overlapping call.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final int SAMPLE_COLUMN = Integer.parseInt(args[1]);  // One-based
        final String OUTPUT_INTERVALS = args[2];
        final String OUTPUT_WINDOWS = args[3];
        final String OUTPUT_GRAPHS = args[4];
        
        boolean isSV;
        char c, d;
        int i;
        int lastInterval, row, pos, length, nCollisions;
        int currentSV, currentChromosome, currentStart, currentEnd, currentGT, currentFirst, sv, chr, start, end, gt;
        String str, tmpString, svtype;
        BufferedReader br;
        BufferedWriter bwIntervals, bwWindows, bwGraphs;
        int[] intervals;
        String[] tokens;
        
        // Loading all VCF records in memory (including records with no GT).
        intervals = new int[5000];  // Arbitrary, multiple of 5.
        lastInterval=-1;
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==VCFconstants.COMMENT) {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            svtype=VCFconstants.getField(tokens[7],VCFconstants.SVTYPE_STR);
            if (svtype!=null && svtype.length()!=0) {
                isSV=true;
                row=svType2Row(svtype);
            }
            else {
                isSV=false;
                row=refAlt2Row(tokens[3],tokens[4]);
            }
            if (row==-1) {
				str=br.readLine();
				continue;
            }
            pos=Integer.parseInt(tokens[1]);
            tmpString=VCFconstants.getField(tokens[7],VCFconstants.SVLEN_STR);
            if (tmpString!=null) length=Integer.parseInt(tmpString);
            else if (row==REPLACEMENT) length=tokens[3].length()-1;
            else length=Math.max(tokens[3].length(),tokens[4].length())-1;
            pos=Integer.parseInt(tokens[1]);
            if (row==DEL || row==INV || row==DUP || row==REPLACEMENT) { start=pos+1; end=pos+length; }
            else if (row==INS) { start=pos; end=pos+1; }
            else if (row==SNP) { start=pos; end=pos; }
            else {
				str=br.readLine();
				continue;
            }
            if (lastInterval+5>=intervals.length) {
                int[] newArray = new int[intervals.length<<1];
                System.arraycopy(intervals,0,newArray,0,intervals.length);
                intervals=newArray;
            }
            intervals[++lastInterval]=isSV?1:0;
            intervals[++lastInterval]=VCFconstants.string2contig(tokens[0]);
            intervals[++lastInterval]=start;  // One-based
            intervals[++lastInterval]=end;  // One-based
            intervals[++lastInterval]=gt2Id(tokens[8+SAMPLE_COLUMN]);
            str=br.readLine();
        }
        br.close();
        if (lastInterval<=4) {
            System.err.println("The VCF contains less than two calls.");
            return;
        }
        bwIntervals = new BufferedWriter(new FileWriter(OUTPUT_INTERVALS));
        for (i=0; i<lastInterval; i+=5) bwIntervals.write(intervals[i]+","+intervals[i+1]+","+intervals[i+2]+","+intervals[i+3]+","+intervals[i+4]+"\n");
        bwIntervals.close();
        
        // Finding windows of overlapping records (including records with 
        // no GT).
        bwWindows = new BufferedWriter(new FileWriter(OUTPUT_WINDOWS));
        bwGraphs = new BufferedWriter(new FileWriter(OUTPUT_GRAPHS));
        bwGraphs.write("graph G {\n");
        currentSV=intervals[0];
        currentChromosome=intervals[1];
        currentStart=intervals[2];
        currentEnd=intervals[3];
        currentGT=intervals[4];
        currentFirst=0;
        for (i=5; i<=lastInterval; i+=5) {
            sv=intervals[i];
            chr=intervals[i+1];
            start=intervals[i+2];
            end=intervals[i+3];
            gt=intervals[i+4];
            if (chr!=currentChromosome || start>currentEnd) {
                buildGraph(intervals,currentFirst,i-1);
                printGraph(intervals,currentFirst,i-1,bwGraphs);
                nCollisions=keepCollisions(intervals,currentFirst,i-1);
                bwWindows.write((currentFirst/5)+","+((i-5)/5)+","+nCollisions+"\n");
                // Next window
                currentSV=sv;
                currentChromosome=chr;
                currentStart=start;
                currentEnd=end;
                currentGT=gt;
                currentFirst=i;
            }
            else if (end>currentEnd) currentEnd=end;
        }
        buildGraph(intervals,currentFirst,lastInterval);
        printGraph(intervals,currentFirst,lastInterval,bwGraphs);
        nCollisions=keepCollisions(intervals,currentFirst,lastInterval);
        bwWindows.write((currentFirst/5)+","+((lastInterval+1-5)/5)+","+nCollisions+"\n");
        bwGraphs.write("}");
        bwWindows.close(); bwGraphs.close();
    }
    
    
	private static final int svType2Row(String type) {
		if ( type.equalsIgnoreCase(VCFconstants.DEL_STR) || 
			 type.equalsIgnoreCase(VCFconstants.DEL_ME_STR)
		   ) return DEL;
		else if (type.equalsIgnoreCase(VCFconstants.INV_STR)) return INV;
        else if ( type.equalsIgnoreCase(VCFconstants.DUP_STR) ||
			      type.equalsIgnoreCase(VCFconstants.DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(VCFconstants.DUP_INT_STR) ||
                  type.equalsIgnoreCase(VCFconstants.CNV_STR)
			    ) return DUP;
        else if ( type.equalsIgnoreCase(VCFconstants.INS_STR) ||
                  type.equalsIgnoreCase(VCFconstants.INS_ME_STR) ||
                  type.equalsIgnoreCase(VCFconstants.INS_NOVEL_STR)
                ) return INS;
		else return -1;
	}
    
    
	private static final int refAlt2Row(String ref, String alt) {
        if (ref.length()==1) {
            if (alt.length()>1) return INS;
            else return SNP;
        }
        else {
            if (alt.length()==1) return DEL;
            else return REPLACEMENT;
        }
	}
    
    
    /**
     * @return a unique ID for every possible genotype.
     */
    private static final int gt2Id(String gt) {
        final char a = gt.charAt(0);
        final char b = gt.charAt(1);
        final char c = gt.charAt(2);
        
        if (b=='/') {
            if (a=='.' || a=='0') {
                if (c=='.' || c=='0') return 0;
                else return 1;
            }
            else {
                if (c=='.' || c=='0') return 2;
                else return 3;
            }
        }
        else {
            if (a=='.' || a=='0') {
                if (c=='.' || c=='0') return 4;
                else return 5;
            }
            else {
                if (c=='.' || c=='0') return 6;
                else return 7;
            }
        }
    }
    
    
    private static final boolean isPresent(int gtId) {
        return gtId!=0 && gtId!=4;
    }
    
    
    private static final boolean isPhased(int gtId) {
        return gtId>=4;
    }
    
    
    private static final String gtId2color(int gtId) {
        switch (gtId) {
            case 0: return "white";
            case 1: return "gray";
            case 2: return "gray";
            case 3: return "gray";
            case 4: return "white";
            case 5: return "blue";
            case 6: return "red";
            case 7: return "violet";
        }
        return "";
    }
    
    
    private static final boolean isCollision(int gtId1, int gtId2) {
        return (gtId1==5 && (gtId2==5 || gtId2==7)) ||
               (gtId1==6 && (gtId2==6 || gtId2==7)) ||
               (gtId1==7 && (gtId2==5 || gtId2==6 || gtId2==7));
    }
    
    
    /**
     * Builds the overlap graph of $intervals[first..last]$, which contains one 
     * node per VCF record, and which contains an edge between two nodes iff the
     * corresponding intervals overlap.
     *
     * @param first,last a range in $intervals$ (inclusive).
     */
    private static final void buildGraph(int[] intervals, int first, int last) throws IOException {
        final int CAPACITY = 10;  // Arbitrary
        
        int i, j;
        int idI, idJ, startI, startJ, endI, endJ;     
        
        // Allocating memory
        nNodes=(last-first+1)/5;
        if (lastNeighbor==null) {
            lastNeighbor = new int[nNodes];
            neighbors = new int[nNodes][CAPACITY];
        }
        else if (lastNeighbor.length<nNodes) {
            lastNeighbor = new int[nNodes];
            int[][] newArray = new int[nNodes][0];
            System.arraycopy(neighbors,0,newArray,0,neighbors.length);
            for (i=neighbors.length; i<nNodes; i++) newArray[i] = new int[CAPACITY];
            neighbors=newArray;
        }
        for (i=0; i<nNodes; i++) lastNeighbor[i]=-1;
      
        // Adding edges
        for (i=first; i<last; i+=5) {
            if (!isPresent(intervals[i+4])) continue;
            idI=(i-first)/5;
            startI=intervals[i+2]; endI=intervals[i+3];
            for (j=i+5; j<last; j+=5) {
                if (!isPresent(intervals[j+4])) continue;
                startJ=intervals[j+2]; endJ=intervals[j+3];
                idJ=(j-first)/5;
                if (startJ>endI) break;
                addNeighbor(idI,idJ);
                addNeighbor(idJ,idI);
            }
        }
        
        // Sorting neighbors
        nEdges=0;
        for (i=0; i<nNodes; i++) {
            Arrays.sort(neighbors[i],0,lastNeighbor[i]+1);
            nEdges+=lastNeighbor[i]+1;
        }
        nEdges>>=1;
    }
    
    
    /**
     * Adds $idJ$ to $neighbors[idI]$.
     */
    private static final void addNeighbor(int idI, int idJ) {
        lastNeighbor[idI]++;
        if (lastNeighbor[idI]==neighbors[idI].length) {
            int[] newArray = new int[neighbors[idI].length<<1];
            System.arraycopy(neighbors[idI],0,newArray,0,neighbors[idI].length);
            neighbors[idI]=newArray;
        }
        neighbors[idI][lastNeighbor[idI]]=idJ;
    }
    
    
    /**
     * Removes every non-collision egde from the graph built by $buildGraph()$.
     *
     * @param first,last a range in $intervals$ (inclusive);
     * @return the number of collisions.
     */
    private static final int keepCollisions(int[] intervals, int first, int last) {
        int i, j, k;
        int neighbor, out;
        
        out=0;
        for (i=0; i<nNodes; i++) {
            k=-1;
            for (j=0; j<=lastNeighbor[i]; j++) {
                neighbor=neighbors[i][j];
                if (isCollision(intervals[first+i*5+4],intervals[first+neighbor*5+4])) neighbors[i][++k]=neighbors[i][j];
            }
            lastNeighbor[i]=k;
            out+=k+1;
        }
        return out>>1;
    }
    
    
    /**
     * Appends to $bw$ a representation of the graph built by $buildGraph()$.
     * Node shape: circle=SNP, box=SV.
     * Node color: genotype.
     * Node label: first position (one-based).
     * Edge: bold=collision, dotted=overlap without collision.
     *
     * Remark: the procedure works both for graphs that contain collision and
     * non-collision edges, and for graphs that contain only collision edges.
     *
     * @param first,last a range in $intervals$ (inclusive).
     */
    private static final void printGraph(int[] intervals, int first, int last, BufferedWriter bw) throws IOException {
        int i, j;
        int neighbor;
        String shape, color;
        
        for (i=0; i<nNodes; i++) {
            if (lastNeighbor[i]==-1) continue;
            shape=intervals[first+i*5]==1?"box":"circle";
            color=gtId2color(intervals[first+i*5+4]);
            bw.write((first/5+i)+" [shape="+shape+",style=filled,fillcolor="+color+",label=\""+intervals[first+i*5+2]+"\"];\n");
        }
        for (i=0; i<nNodes; i++) {   
            for (j=0; j<=lastNeighbor[i]; j++) {
                neighbor=neighbors[i][j];
                if (neighbor>i) {
                    if (isCollision(intervals[first+i*5+4],intervals[first+neighbor*5+4])) bw.write((first/5+i)+" -- "+(first/5+neighbor)+" [style=bold];\n");
                    else bw.write((first/5+i)+" -- "+(first/5+neighbor)+" [style=dotted];\n");
                }
            }
        }
    }
    
    
    /**
     * Sorts nodes by increasing last position
     */
    private static final void sortNodes(int[] intervals, int first, int last) {
        int i;
        int endI, endJ, nodeID;
        
        if (sorted==null || sorted.length==0) {
            sorted = new Interval[nNodes];
            for (i=0; i<nNodes; i++) sorted[i] = new Interval();
        }
        else if (nNodes>sorted.length) {
            Interval[] newArray = new Interval[nNodes];
            System.arraycopy(sorted,0,newArray,0,sorted.length);
            for (i=sorted.length; i<nNodes; i++) newArray[i] = new Interval();
            sorted=newArray;
        }
        for (i=0; i<nNodes; i++) {
            sorted[i].id=i;
            sorted[i].first=intervals[first+i*5+2];
            sorted[i].last=intervals[first+i*5+3];
        }
        Arrays.sort(sorted,0,nNodes);
        if (node2sorted==null || node2sorted.length<nNodes) node2sorted = new int[nNodes];
        for (i=0; i<nNodes; i++) node2sorted[sorted[i].id]=i;
    }
    
    
    /**
     * Computes a maximum-size independent set in the collision graph (using a
     * greedy algorithm that is optimal for interval graphs).
     *
     * @param out output array, indicates which nodes of the graph to keep.
     */
    private static final void maxSizeIndependentSet(int[] intervals, int first, int last, boolean[] out) {
        int i, j;
        int endI, endJ, nodeID;
        
        sortNodes(intervals,first,last);
        Arrays.fill(out,0,nNodes,true);
        for (i=0; i<nNodes; i++) {
            nodeID=sorted[i].id;
            if (!out[nodeID]) continue;
            for (j=0; j<=lastNeighbor[nodeID]; j++) out[neighbors[nodeID][j]]=false;
        }
    }
    
    
    
    /**
     * Computes a maximum-weight independent set in the collision graph.
     *
     *
     * @param weights one positive score per node.
     */
    private static final void maxScoreIndependentSet(int[] intervals, int first, int last, double[] weights, boolean[] out) {
        
        
        
        sortNodes(intervals,first,last);
        Arrays.fill(out,0,nNodes,true);
        
        
        
    }
    
    
    
    
    private static class Interval implements Comparable {
        public int id, first, last;
        
        public boolean equals(Object other) {
            final Interval otherInterval = (Interval)other;
            return last==otherInterval.last;
        }
        
        public int compareTo(Object other) {
            final Interval otherInterval = (Interval)other;
            if (last<otherInterval.last) return -1;
            else if (last>otherInterval.last) return 1;
            else return 0;
        }
    }
    
    
    

}