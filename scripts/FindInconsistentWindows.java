import java.util.Arrays;
import java.io.*;


/**
 * 
 */
public class FindInconsistentWindows {
    /**
     * Temporary variables to store the graph of a window
     */
    private static int nNodes, nEdges;
    private static int[] lastNeighbor, interval2id, color, queue;
    private static int[][] neighbors;
    
    
    /**
     * @param args 0 the VCF is assumed to be sorted.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final String OUTPUT_INTERVALS = args[1];
        final String OUTPUT_BIPARTITE = args[2];
        final String OUTPUT_NOT_BIPARTITE = args[3];
        
        int i;
        int lastInterval, row, pos, length;
        int currentChromosome, currentStart, currentEnd, currentFirst, chr, start, end;
        String str, tmpString;
        BufferedReader br;
        BufferedWriter bwIntervals, bwBipartite, bwNotBipartite;
        int[] intervals;
        String[] tokens;
        
        // Loading all intervals, assuming that the VCF is sorted.
        intervals = new int[4000];  // Arbitrary, multiple of 4.
        lastInterval=-1;
        br = new BufferedReader(new FileReader(INPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==VCFconstants.COMMENT) {
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
            pos=Integer.parseInt(tokens[1]);
            if (row==0 || row==1 || row==2) end=pos+length;
            else if (row==3) end=pos+1;
            else {
				str=br.readLine();
				continue;
            }
            if (lastInterval+4>=intervals.length) {
                int[] newArray = new int[intervals.length<<1];
                System.arraycopy(intervals,0,newArray,0,intervals.length);
                intervals=newArray;
            }
            intervals[++lastInterval]=VCFconstants.string2contig(tokens[0]);
            intervals[++lastInterval]=pos+1;  // One-based
            intervals[++lastInterval]=end;
            intervals[++lastInterval]=(tokens[9].charAt(0)=='1'?1:0)+(tokens[9].charAt(2)=='1'?1:0);
            str=br.readLine();
        }
        br.close();
        if (lastInterval<=3) {
            System.err.println("The VCF contains just one call.");
            return;
        }
        bwIntervals = new BufferedWriter(new FileWriter(OUTPUT_INTERVALS));
        for (i=0; i<lastInterval; i+=4) bwIntervals.write(intervals[i]+","+intervals[i+1]+","+intervals[i+2]+","+intervals[i+3]+"\n");
        bwIntervals.close();
        
        // Finding windows of overlapping intervals
        bwBipartite = new BufferedWriter(new FileWriter(OUTPUT_BIPARTITE));
        bwNotBipartite = new BufferedWriter(new FileWriter(OUTPUT_NOT_BIPARTITE));
        currentChromosome=intervals[0];
        currentStart=intervals[1];
        currentEnd=intervals[2];
        currentFirst=0;
        for (i=4; i<=lastInterval; i+=4) {
            chr=intervals[i];
            start=intervals[i+1];
            end=intervals[i+2];
            if (chr!=currentChromosome || start>currentEnd) {
                if (currentFirst==i-4) bwBipartite.write("1,"+intervals[i-1]+","+(currentFirst>>2)+","+(currentFirst>>2)+"\n");
                else {
                    buildGraph(intervals,currentFirst,i-1);
                    if (isBipartite()) bwBipartite.write(((i-currentFirst)>>2)+","+nNodes+","+(currentFirst>>2)+","+((i-4)>>2)+"\n");
                    else bwNotBipartite.write(((i-currentFirst)>>2)+","+nNodes+","+(currentFirst>>2)+","+((i-4)>>2)+"\n");
                }
                // Next window
                currentChromosome=chr;
                currentStart=start;
                currentEnd=end;
                currentFirst=i;
            }
            else if (end>currentEnd) currentEnd=end;
        }
        if (currentFirst==lastInterval-3) bwBipartite.write("1,"+intervals[lastInterval]+","+(currentFirst>>2)+","+(currentFirst>>2)+"\n");
        else {
            buildGraph(intervals,currentFirst,lastInterval);
            if (isBipartite()) bwBipartite.write(((lastInterval+1-currentFirst)>>2)+","+nNodes+","+(currentFirst>>2)+","+((lastInterval-3)>>2)+"\n");
            else bwNotBipartite.write(((lastInterval+1-currentFirst)>>2)+","+nNodes+","+(currentFirst>>2)+","+((lastInterval-3)>>2)+"\n");
        }
        bwBipartite.close(); bwNotBipartite.close();
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
    
    
    /**
     * Builds the overlap graph of $intervals[first..last]$, which contains a 
     * vertex for each instance of a call that appears in its GT, and an edge 
     * between two nodes iff the corresponding intervals intersect.
     */
    private static final void buildGraph(int[] intervals, int first, int last) throws IOException {
        final int CAPACITY = 10;  // Arbitrary
        int i, j;
        int idI, idJ, startI, startJ, endI, endJ;        
        
        // Allocating memory
        if (interval2id==null || interval2id.length<(last-first+1)>>2) interval2id = new int[(last-first+1)>>2];
        nNodes=0;
        for (i=first; i<=last; i+=4) {
            interval2id[(i-first)>>2]=nNodes;
            nNodes+=intervals[i+3];
        }
        if (lastNeighbor==null) {
            lastNeighbor = new int[nNodes];
            for (i=0; i<nNodes; i++) lastNeighbor[i]=-1;
            neighbors = new int[nNodes][CAPACITY];
        }
        else if (lastNeighbor.length<nNodes) {
            lastNeighbor = new int[nNodes];
            for (i=0; i<nNodes; i++) lastNeighbor[i]=-1;
            int[][] newArray = new int[nNodes][0];
            System.arraycopy(neighbors,0,newArray,0,neighbors.length);
            for (i=neighbors.length; i<nNodes; i++) newArray[i] = new int[CAPACITY];
            neighbors=newArray;
        }
        
        // Adding self-edges, if any.
        for (i=first; i<last; i+=4) {
            if (intervals[i+3]!=2) continue;
            idI=interval2id[(i-first)>>2];
            addNeighbor(idI,idI+1);
            addNeighbor(idI+1,idI);
        }
        
        // Adding non-self edges, if any.
        for (i=first; i<last; i+=4) {
            if (intervals[i+3]==0) continue;
            idI=interval2id[(i-first)>>2];
            startI=intervals[i+1]; endI=intervals[i+2];
            for (j=i+4; j<last; j+=4) {
                if (intervals[j+3]==0) continue;
                startJ=intervals[j+1]; endJ=intervals[j+2];
                idJ=interval2id[(j-first)>>2];
                if (startJ>endI) break;
                if (intervals[i+3]==1) {
                    if (intervals[j+3]==1) {
                        addNeighbor(idI,idJ);
                        addNeighbor(idJ,idI);
                    }
                    else {
                        addNeighbor(idI,idJ); addNeighbor(idI,idJ+1);
                        addNeighbor(idJ,idI); addNeighbor(idJ+1,idI);
                    }
                }
                else {
                    if (intervals[j+3]==1) {
                        addNeighbor(idI,idJ); addNeighbor(idI+1,idJ);
                        addNeighbor(idJ,idI); addNeighbor(idJ,idI+1);
                    }
                    else {
                        addNeighbor(idI,idJ); addNeighbor(idI+1,idJ); 
                        addNeighbor(idI,idJ+1); addNeighbor(idI+1,idJ+1);
                        addNeighbor(idJ,idI); addNeighbor(idJ,idI+1);
                        addNeighbor(idJ+1,idI); addNeighbor(idJ+1,idI+1);
                    }
                }
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
     * Simple breadth-first search test.
     * Works on the graph currently stored in the global variables.
     */
    private static final boolean isBipartite() {
        int i, j;
        int first, last, currentNode, currentColor, neighbor;
        
        if (nEdges<=2) return true;
        if (color==null || color.length<nNodes) color = new int[nNodes];
        for (i=0; i<nNodes; i++) color[i]=-1;
        if (queue==null || queue.length<nNodes<<1) queue = new int[nNodes<<1];
        for (i=0; i<nNodes; i++) {
            if (color[i]!=-1) continue;
            color[i]=0; queue[0]=i; first=0; last=0; 
            while (first<=last) {
                currentNode=queue[first];
                first=(first+1)%queue.length;
                currentColor=color[currentNode];
                for (j=0; j<=lastNeighbor[currentNode]; j++) {
                    neighbor=neighbors[currentNode][j];
                    if (color[neighbor]==-1) {
                        color[neighbor]=currentColor==0?1:0;
                        last=(last+1)%queue.length;
                        queue[last]=neighbor;
                    }
                    else if (color[neighbor]==currentColor) return false;
                }
            }
        }
        return true;
    }

}