import java.util.Arrays;
import java.util.Random;
import java.io.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.imageio.*;
import java.text.*;


/**
 * 
 */
public class DrawPhasedVariants {
    /**
     * Expected weight tags in the INFO field and in the SAMPLE field of a VCF
     */
    private static final String WEIGHT_ID_INFO = "IS_WEIGHT";
    private static final String WEIGHT_ID_SAMPLE = "IS_WEIGHT";
    private static boolean LOAD_WEIGHT_FROM_SAMPLE_COLUMN;
    
    /**
     * Variant types
     */
    private static final int DEL = 0;
    private static final int INV = 1;
    private static final int DUP = 2;
    private static final int INS = 3;
    private static final int SNP = 4;
    private static final int REPLACEMENT = 5;
    
    /**
     * Genotypes
     */
    private static final int PHASED_00 = 0;
    private static final int PHASED_01 = 1;
    private static final int PHASED_10 = 2;
    private static final int PHASED_11 = 3;
    private static final int UNPHASED_00 = 4;
    private static final int UNPHASED_01 = 5;
    private static final int UNPHASED_10 = 6;
    private static final int UNPHASED_11 = 7;
    
    /**
     * Drawing constants
     */
    private static final int PIXELS_PER_POS = 50;
    private static final int PIXELS_DELTA = 5;
    private static final int N_ROWS = 5*PIXELS_PER_POS;
    private static final Color COLOR_BACKGROUND = new Color(0x00FFFFFF);
    private static final Stroke STROKE_UNPHASED = new BasicStroke(1.0f,BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND,10.0f,new float[] {10.0f},0.0f);
    private static final Stroke STROKE_PHASED = new BasicStroke(1.0f,BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND,10.0f);
    private static final Stroke STROKE_BACKGROUND = new BasicStroke(0.0f);
    private static Color[] type2color;
    
    /**
     * A maximal set of overlapping intervals.
     *
     * Remark: the last interval in $window$ might not overlap with the previous
     * intervals, but it might be the first interval of the next set of 
     * overlapping intervals.
     */
    private static Interval[] window;
    
    /**
     * Last element in $window$.
     */
    private static int windowLast;
    
    /**
     * Last position of an overlapping interval in $window$ (one-based,
     * inclusive). Might be smaller than the last position of the last interval
     * in $window$, if such an interval does not overlap with the previous ones.
     */
    private static int windowLastPos;
    
    

    /**
     * @param args
     * 2: Only build the collision histogram, do not draw (1/0).
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final String OUTPUT_DIR = args[1];
        final boolean COUNT_ONLY = Integer.parseInt(args[2])==1;
        LOAD_WEIGHT_FROM_SAMPLE_COLUMN = Integer.parseInt(args[3])==1;
        
        int i;
        BufferedReader br;
        long[] histogram;
        
        initType2color();
        window = new Interval[100];
        for (i=0; i<window.length; i++) window[i] = new Interval();
        histogram = new long[100];
        windowLast=-1; windowLastPos=-1;
        br = new BufferedReader(new FileReader(INPUT_VCF));
        while (true) {
            initNextWindow();
System.err.println("VITTU> 1  windowLast="+windowLast+" windowLastPos="+windowLastPos);            
            loadNextWindow(br);
System.err.println("VITTU> 2  windowLast="+windowLast+" windowLastPos="+windowLastPos+(windowLast>=0?" first pos: "+window[0].first:""));
            if (windowLast==-1) break;
            drawWindow(OUTPUT_DIR,histogram,COUNT_ONLY);
        }
        br.close();
        System.err.println("Histogram (nCollision, nHaplotypes):");
        for (i=0; i<histogram.length; i++) System.err.println(i+","+histogram[i]);
    }
    
    
    
    
    
    
    
    
    // --------------------- INDEPENDENT SET PROCEDURES ------------------------
    
    /**
     * The optimal solution of independent set procedures
     */
    private static double maxWeight;
    private static boolean[] maxPath;  // <-------- MIGHT NOT BE NEEDED, CAN REUSE THE BOOLEAN INSIDE INTERVAL?
    
    /**
     * Space reused by independent set procedures
     */
    private static Interval[] sampleWindow;
    private static int sampleWindowLast;
    private static Point[] endpoints;
    private static int lastEndpoint;
    
    
    /**
     * Assume that we can fix the genotypes of $sample$ only by setting some
     * of them to 0/0, i.e. by deleting records on both haplotypes. Then we
     * would like to keep a heaviest set of VCF records that do not overlap on
     * any haplotype. The procedure computes a max-weight independent set of the 
     * graph where every element of $window[0..last]$ that occurs in $sample$ 
     * (on any haplotype) is a node, and where there is an edge iff two elements
     * overlap on some haplotype.
     *
     * Remark: the procedure considers also unphased calls, since they might
     * collide with phased calls.
     *
     * Remark: the algorithm is a simple tree search on the (exponential) space
     * of all possible solutions, with a simple pruning based on the current
     * best solution.
     * 
     * Remark: the procedure assumes that $window[0..last]$ is sorted by last
     * position.
     */
    private static final void independentSet1(int sample, int last) {
        boolean error;
        int i, j;
        double weight, currentWeight;
        Interval currentInterval, neighbor;
        boolean[] active;

        // Collecting only the intervals that are present in $sample$ (some of
        // which might be unphased), and assigning a weight to them.
        if (sampleWindow==null || sampleWindow.length<last+1) sampleWindow = new Interval[last+1];
        sampleWindowLast=-1;
        for (i=0; i<=last; i++) {
            currentInterval=window[i];
            if (!currentInterval.isPresent(sample)) continue;
            sampleWindow[++sampleWindowLast]=currentInterval;
            currentInterval.setWeight(sample);
        }
        
        // Recursion
        maxWeight=0;
        if (maxPath==null || maxPath.length<sampleWindowLast+1) maxPath = new boolean[sampleWindowLast+1];
        Arrays.fill(maxPath,false);
        active = new boolean[sampleWindowLast+1];
        Arrays.fill(active,true);
        for (i=sampleWindowLast; i>=0; i--) independendSet1_impl(i,active,new int[0],0,sample);
        error=false;
        if (maxWeight==0) error=true;
        else {
            error=true;
            for (i=0; i<maxPath.length; i++) {
                if (maxPath[i]) { error=false; break; }
            }
        }
        if (error) {
            System.err.println("No selected interval in the optimal solution of sample "+sample+"?!  maxWeight="+maxWeight+" sampleWindowLast="+sampleWindowLast);        
            System.exit(1);
        }
        
        // Removing from $sample$ every phased interval that is not in the
        // max-weight independent set.
        for (i=0; i<=sampleWindowLast; i++) {
            if (!maxPath[i]) sampleWindow[i].genotypes[sample]=UNPHASED_00;
        }
    }
    
    
    /**
     * Remark: the procedure assumes that $maxPath$ has already been allocated. 
     *
     * @param id position in $sampleWindow$ to solve for;
     * @param active one flag per position in $sampleWindow$; indicates whether 
     * the corresponding interval is active before solving for $id$;
     * @param ancestors ancestors of $id$ in the search tree;
     * @param ancestorsWeight sum of weights of all the elements in $ancestors$.
     */
    private static final void independendSet1_impl(int id, boolean[] active, int[] ancestors, double ancestorsWeight, int sample) {
        final Interval currentInterval = sampleWindow[id];
        final double weight_prime = ancestorsWeight+currentInterval.weight;
        
        boolean hasChild;
        int i;
        double upperBound;
        boolean[] active_prime;
        int[] ancestors_prime;
        
        active_prime = new boolean[active.length];
        System.arraycopy(active,0,active_prime,0,active.length);
        active_prime[id]=false;
        hasChild=false; upperBound=weight_prime;
        for (i=id-1; i>=0; i--) {
            if (active_prime[i] && !sampleWindow[i].precedes(currentInterval,sample)) active_prime[i]=false;
            if (active_prime[i]) {
                hasChild=true;
                upperBound+=sampleWindow[i].weight;
            }
        }
        if (upperBound<maxWeight) return;  // Pruning all recursive calls
        if (hasChild) {
            ancestors_prime = new int[ancestors.length+1];
            System.arraycopy(ancestors,0,ancestors_prime,0,ancestors.length);
            ancestors_prime[ancestors_prime.length-1]=id;
            for (i=id-1; i>=0; i--) {
                if (!active_prime[i]) continue;
                if (upperBound<maxWeight) break;  // Pruning this and all the following recursive calls
                independendSet1_impl(i,active_prime,ancestors_prime,weight_prime,sample);
                upperBound-=sampleWindow[i].weight;
            }
        }
        else if (weight_prime>maxWeight) {  // Base case of the recursion
            maxWeight=weight_prime;
            Arrays.fill(maxPath,false);
            for (i=0; i<ancestors.length; i++) maxPath[ancestors[i]]=true;
            maxPath[id]=true;
        }
    }
    
    
    /**
     * Assume that we can fix the genotypes of $sample$ only by replacing ones 
     * with zeros, i.e. by deleting specific calls on specific haplotypes. Then 
     * we would like to keep, on each haplotype, a heaviest set of calls that do
     * not overlap on that haplotype. The procedure computes, for each 
     * haplotype, a max-weight independent set of the graph where every element 
     * of $window[0..last]$ that occurs on that haplotype is a node, and where 
     * there is an edge iff two elements overlap.
     *
     * Remark: the procedure considers only phased calls and 1/1 unphased calls.
     *
     * Remark: the procedure assumes that $window[0..last]$ is sorted by last
     * position.
     */
    private static final void independentSet2(int sample, int last) {
        int i;
        
        // Allocating reused space
        if (sampleWindow==null || sampleWindow.length<last+1) sampleWindow = new Interval[last+1];
        if (endpoints==null) {
            endpoints = new Point[(last+1)<<1];
            for (i=0; i<endpoints.length; i++) endpoints[i] = new Point();
        }
        else if (endpoints.length<(last+1)<<1) {
            Point[] newArray = new Point[(last+1)<<1];
            System.arraycopy(endpoints,0,newArray,0,endpoints.length);
            for (i=endpoints.length; i<newArray.length; i++) newArray[i] = new Point();
            endpoints=newArray;
        }
        
        // Solving Hap 1
        sampleWindowLast=-1;
        for (i=0; i<=last; i++) {
            if (window[i].onHap1(sample)) {
                sampleWindow[++sampleWindowLast]=window[i];
                window[i].setWeight(sample);
            }
        }
        independentSet2_impl(sample,true);
        
        // Solving Hap 2
        sampleWindowLast=-1;
        for (i=0; i<=last; i++) {
            if (window[i].onHap2(sample)) {
                sampleWindow[++sampleWindowLast]=window[i];
                window[i].setWeight(sample);
            }
        }
        independentSet2_impl(sample,false);
    }


    /**
     * Implements the linear algorithm in Section 2 of:
     *
     * Hsiao et al. "An efficient algorithm for finding a maximum weight
     * 2-independent set on interval graphs." Information Processing Letters
     * 43.5 (1992): 229-235.
     *
     * @param onHap1 the procedure assumes that all intervals in $sampleWindow$ 
     * occur on haplotype 1 (TRUE) or on haplotype 2 (FALSE) of $sample$.
     */
    private static final void independentSet2_impl(int sample, boolean onHap1) {
        boolean error;
        int i, j;
        double availableWeight;
        Interval currentInterval, previous;

        // Computing a max-weight independent set
        intervals2endpoints();
        for (i=0; i<=sampleWindowLast; i++) sampleWindow[i].clearIndependentSetVariables();
        availableWeight=0.0; previous=null;
        for (i=0; i<=lastEndpoint; i++) {
            for (j=0; j<=endpoints[i].lastOpen; j++) {
                currentInterval=endpoints[i].open[j];
                currentInterval.independentSetWeight=availableWeight+currentInterval.weight;
                currentInterval.independentSetPrevious=previous;
            }
            for (j=0; j<=endpoints[i].lastClosed; j++) {
                currentInterval=endpoints[i].closed[j];
                if (currentInterval.independentSetWeight>availableWeight) {
                    availableWeight=currentInterval.independentSetWeight;
                    previous=currentInterval;
                }
            }
        }
        maxWeight=availableWeight;
        if (maxWeight==0) {
            System.err.println("The optimal solution of sample "+sample+" "+(onHap1?"hap1":"hap2")+" has weight zero?! sampleWindowLast="+sampleWindowLast);        
            System.exit(1);
        }

        // Removing from the current haplotype every interval that is not in the
        // selected max-weight independent set.
        for (i=sampleWindowLast; i>=0; i--) {
            if (sampleWindow[i].independentSetWeight!=maxWeight) continue;
            currentInterval=sampleWindow[i];
            while (currentInterval!=null) {
                currentInterval.inIndependentSet=true;
                currentInterval=currentInterval.independentSetPrevious;
            }
            break;
        }
        error=true;
        if (onHap1) {
            for (i=0; i<=sampleWindowLast; i++) {
                if (!sampleWindow[i].inIndependentSet) sampleWindow[i].genotypes[sample]=REMOVE_HAP1[sampleWindow[i].genotypes[sample]];
                else error=false;
            }
        }
        else {
            for (i=0; i<=sampleWindowLast; i++) {
                if (!sampleWindow[i].inIndependentSet) sampleWindow[i].genotypes[sample]=REMOVE_HAP2[sampleWindow[i].genotypes[sample]];
                else error=false;
            }
        }
        if (error) {
            System.err.println("No selected interval in the optimal solution of sample "+sample+" "+(onHap1?"hap1":"hap2")+"?!  maxWeight="+maxWeight+" sampleWindowLast="+sampleWindowLast);        
            System.exit(1);
        }
    }
    
    
    private static final int[] REMOVE_HAP1 = new int[] {
        // 0|0, 0|1, 1|0, 1|1,  0/0, 0/1, 1/0, 1/1
        PHASED_00,PHASED_01,PHASED_00,PHASED_01,  UNPHASED_00,UNPHASED_01,UNPHASED_10,UNPHASED_01
    };
    
    
    private static final int[] REMOVE_HAP2 = new int[] {
        // 0|0, 0|1, 1|0, 1|1,  0/0, 0/1, 1/0, 1/1
        PHASED_00,PHASED_00,PHASED_10,PHASED_10,  UNPHASED_00,UNPHASED_01,UNPHASED_10,UNPHASED_01
    };


    /**
     * Stores in $endpoints$ the sorted list of distinct first and last 
     * positions of all intervals in $sampleWindow$.
     *
     * Remark: the procedure assumes that $endpoints$ has already been 
     * initialized.
     */
    private static final void intervals2endpoints() {
        int i, j;
        Point tmpPoint;

        lastEndpoint=-1;
        for (i=0; i<=sampleWindowLast; i++) {
            lastEndpoint++;
            endpoints[lastEndpoint].clear();
            endpoints[lastEndpoint].pos=sampleWindow[i].first;
            endpoints[lastEndpoint].addOpen(sampleWindow[i]);
            lastEndpoint++;
            endpoints[lastEndpoint].clear();
            endpoints[lastEndpoint].pos=sampleWindow[i].last;
            endpoints[lastEndpoint].addClosed(sampleWindow[i]);
        }
        if (lastEndpoint==1) return;
        Arrays.sort(endpoints,0,lastEndpoint+1);
        j=0;
        for (i=1; i<=lastEndpoint; i++) {
            if (endpoints[i].pos==endpoints[j].pos) endpoints[j].merge(endpoints[i]);
            else {
                j++;
                tmpPoint=endpoints[j];
                endpoints[j]=endpoints[i];
                endpoints[i]=tmpPoint;
            }
        }
        lastEndpoint=j;
    }


    /**
     * The endpoint of an interval
     */
    private static class Point implements Comparable {
        private static final int CAPACITY = 2;  // Arbitrary
        public int pos;  // One-based
        public Interval[] open, closed;
        public int lastOpen, lastClosed;


        public Point() {
            open = new Interval[CAPACITY];
            closed = new Interval[CAPACITY];
            clear();
        }


        public final void clear() {
            pos=-1; lastOpen=-1; lastClosed=-1;
        }


        public final void addOpen(Interval interval) {
            lastOpen++;
            if (lastOpen==open.length) {
                Interval[] newArray = new Interval[open.length<<1];
                System.arraycopy(open,0,newArray,0,open.length);
                open=newArray;
            }
            open[lastOpen]=interval;
        }


        public final void addClosed(Interval interval) {
            lastClosed++;
            if (lastClosed==closed.length) {
                Interval[] newArray = new Interval[closed.length<<1];
                System.arraycopy(closed,0,newArray,0,closed.length);
                closed=newArray;
            }
            closed[lastClosed]=interval;
        }


        /**
         * Adds to $open,closed$ the intervals in the corresponding arrays of
         * $otherPoint$.
         */
        public final void merge(Point otherPoint) {
            int newLength;

            newLength=lastOpen+1+otherPoint.lastOpen+1;
            if (newLength>open.length) {
                Interval[] newArray = new Interval[newLength];
                System.arraycopy(open,0,newArray,0,lastOpen+1);
                open=newArray;
            }
            System.arraycopy(otherPoint.open,0,open,lastOpen+1,otherPoint.lastOpen+1);
            lastOpen=newLength-1;
            newLength=lastClosed+1+otherPoint.lastClosed+1;
            if (newLength>closed.length) {
                Interval[] newArray = new Interval[newLength];
                System.arraycopy(closed,0,newArray,0,lastClosed+1);
                closed=newArray;
            }
            System.arraycopy(otherPoint.closed,0,closed,lastClosed+1,otherPoint.lastClosed+1);
            lastClosed=newLength-1;
        }


        public boolean equals(Object other) {
            final Point otherPoint = (Point)other;
            return pos==otherPoint.pos;
        }


        public int compareTo(Object other) {
            final Point otherPoint = (Point)other;
            if (pos<otherPoint.pos) return -1;
            else if (pos>otherPoint.pos) return 1;
            return 0;
        }


        public String toString() {
            int i;
            String out;

            out="pos="+pos+"\n  open: ";
            for (i=0; i<=lastOpen; i++) out+=open[i].variantType+",";
            out+="\n  closed: ";
            for (i=0; i<=lastClosed; i++) out+=closed[i].variantType+",";
            return out+"\n";
        }
    }

    
    
    
    // --------------------- WINDOW DRAWING PROCEDURES -------------------------
    
    /**
     * Stores in a directory the overlap diagram of the window for every sample
     * that contains a collision (one PNG file per sample).
     *
     * @param histogram for each $i$, the number of haplotypes (=windows in 
     * samples) with $i$ collisions. Windows with exactly one variant are not
     * considered in the count.
     * @param countOnly only compute the histogram, do not draw.
     */
    private static final void drawWindow(String outputDir, long[] histogram, boolean countOnly) throws IOException {
        final int FIRST = window[0].first;
        final int N_COLUMNS = (2+windowLastPos-FIRST+1)*PIXELS_PER_POS;
        final int LAST = isLastInWindow()?windowLast:windowLast-1;
        final int N_SAMPLES = window[0].genotypes.length;
        final int CHR = window[0].chr;
        final String WINDOW_DIR = outputDir+"/chr"+CHR+"_"+FIRST+"_"+windowLastPos;
        final int HIGH_COLLISIONS = 10;  // Arbitrary
        
        int i, j, x, y;
        int gt, width, collisions;
        Random random = null;
        File directory = null;
        BufferedImage image = null;
        Graphics2D graphics = null;
        
        if (LAST==0) return;  // Ignoring windows with only one variant
        Arrays.sort(window,0,LAST+1);
        if (!countOnly) {
            random = new Random();
            directory = new File(WINDOW_DIR);
            image = new BufferedImage(N_COLUMNS,N_ROWS,BufferedImage.TYPE_INT_RGB);
            graphics=image.createGraphics();
        }
        for (j=0; j<N_SAMPLES; j++) {
            collisions=nCollisions(LAST,j);
            histogram[collisions>histogram.length-1?histogram.length-1:collisions]++;
            if (collisions==0 || countOnly) continue;
            if (collisions>=HIGH_COLLISIONS) System.err.println(collisions+" collisions in window "+WINDOW_DIR+" sample "+j);
            if (!directory.exists()) directory.mkdirs();
            drawWindow_impl(j,LAST,graphics,N_COLUMNS,FIRST,random);
            ImageIO.write(image,"png",new File(WINDOW_DIR+"/sample"+j+"_before.png"));
            independentSet2(j,LAST);
            drawWindow_impl(j,LAST,graphics,N_COLUMNS,FIRST,random);
            ImageIO.write(image,"png",new File(WINDOW_DIR+"/sample"+j+"_after.png"));
        }
    }
    
    
    private static final void drawWindow_impl(int sample, int last, Graphics2D graphics, int nColumns, int first, Random random) {
        int i, x, y;
        int gt, width;
        
        graphics.setStroke(STROKE_BACKGROUND);
        graphics.setColor(COLOR_BACKGROUND);
        graphics.fillRect(0,0,nColumns,N_ROWS);
        for (i=0; i<=last; i++) {
            gt=window[i].genotypes[sample];
            if (gt==UNPHASED_00 || gt==PHASED_00) continue;
            if (gt==UNPHASED_01 || gt==UNPHASED_10) graphics.setStroke(STROKE_UNPHASED);
            else graphics.setStroke(STROKE_PHASED);
            graphics.setColor(type2color[window[i].variantType]);
            x=(1+window[i].first-first)*PIXELS_PER_POS-PIXELS_DELTA+random.nextInt((PIXELS_DELTA)<<1);
            width=(window[i].last-window[i].first+1)*PIXELS_PER_POS;
            if (gt==UNPHASED_01 || gt==UNPHASED_10 || gt==UNPHASED_11 || gt==PHASED_11) {
                y=PIXELS_PER_POS-PIXELS_DELTA+random.nextInt((PIXELS_DELTA)<<1);
                graphics.drawRect(x,y,width,3*PIXELS_PER_POS);
            }
            else if (gt==PHASED_10) {
                y=PIXELS_PER_POS-PIXELS_DELTA+random.nextInt((PIXELS_DELTA)<<1);
                graphics.drawRect(x,y,width,PIXELS_PER_POS);
            }
            else if (gt==PHASED_01) {
                y=PIXELS_PER_POS*3-PIXELS_DELTA+random.nextInt((PIXELS_DELTA)<<1);
                graphics.drawRect(x,y,width,PIXELS_PER_POS);
            }
        }
    }
    
    
    /**
     * Remark: the procedure assumes that $window$ is sorted by last position.
     *
     * @return the number of pairs of intervals in $window[0..last]$ that 
     * collide in $sample$, assuming that a diploid call consists of two 
     * intervals.
     */
    private static final int nCollisions(int last, int sample) {
        int i, j;
        int startI, endI, gtI, out;
        
        out=0;
        for (i=last; i>=0; i--) {
            if (!window[i].isPresent(sample)) continue;
            startI=window[i].first; gtI=window[i].genotypes[sample];
            for (j=i-1; j>=0; j--) {
                if (window[j].last<startI) break;
                out+=N_GT_COLLISIONS[window[j].genotypes[sample]][gtI];
            }
        }
        return out;
    }
    
    
    /**
     * Rows, columns: GT ids. Cells: {0,1,2} the number of pairs of intervals
     * that collide, assuming that a diploid call consists of two intervals. 
     */
    private static final int[][] N_GT_COLLISIONS = new int[][] { 
        // 0|0, 0|1, 1|0, 1|1,  0/0, 0/1, 1/0, 1/1
        {0,0,0,0, 0,0,0,0},
        {0,1,0,1, 0,0,0,1},
        {0,0,1,1, 0,0,0,1},
        {0,1,1,2, 0,0,0,2},
        {0,0,0,0, 0,0,0,0},
        {0,0,0,1, 0,0,0,1},
        {0,0,0,1, 0,0,0,1},
        {0,1,1,2, 0,0,0,2},
    };

    
    
    
    // --------------------- WINDOW LOADING PROCEDURES -------------------------
    
    /**
     * @return TRUE iff the last interval of $window$ overlaps with the previous
     * intervals.
     */
    private static final boolean isLastInWindow() {
        return windowLast==-1 || windowLastPos==-1 || (window[windowLast].chr==window[0].chr && window[windowLast].first<=windowLastPos);
    }
    
    
    /**
     * Prepares $window$ for loading the next set of overlapping intervals.
     */
    private static final void initNextWindow() {
        if (isLastInWindow()) {
            windowLast=-1;
            windowLastPos=-1;
        }
        else {
            Interval tmpInterval = window[0];
            window[0]=window[windowLast];
            window[windowLast]=tmpInterval;
            windowLastPos=window[0].last;
            windowLast=0;
        }
    }
    
    
    /**
     * Loads in $window$ a maximal set of overlapping intervals, regardless of
     * their genotypes.
     *
     * Remark: after this procedure completes, the last interval in $window$
     * might not overlap with the previous ones.
     */
    private static final void loadNextWindow(BufferedReader br) throws IOException {
        int i;
        String str;
        Interval tmpInterval;
        
        while (true) {
            str=br.readLine();
            if (str==null) break;
            if (str.charAt(0)==VCFconstants.COMMENT) continue;
            windowLast++;
            if (windowLast==window.length) {
                Interval[] newArray = new Interval[window.length<<1];
                System.arraycopy(window,0,newArray,0,window.length);
                for (i=window.length; i<newArray.length; i++) newArray[i] = new Interval();
                window=newArray;
            }
            tmpInterval=window[windowLast];
            if (!tmpInterval.init(str)) {
                windowLast--;
                str=br.readLine();
                continue;
            }
            if (!isLastInWindow()) break;
            if (tmpInterval.last>windowLastPos) windowLastPos=tmpInterval.last;
        }
        
        
System.err.println("loadNextWindow> loaded window:");
for (int x=0; x<=windowLast; x++) System.err.println(window[x]);
        
        
        
    }
    
    
    
    
    // ----------------------- PROPERTIES OF INTERVALS  ------------------------
    
    /**
     * A VCF record represented as a weighted interval on the line.
     * The same object is intended to be reused by multiple VCF records.
     */
    private static class Interval implements Comparable {
        private static final String VCF_SEPARATOR = "\t";
        private static final char GT_SEPARATOR = ':';
        
        public boolean isSV;  // True=SV, False=SNP/small indel.
        public int variantType;
        public int chr;
        public int first, last;  // One-based, inclusive.
        public double weight;
        public String[] vcfRecord;  // One cell per VCF column
        public int[] genotypes;
        
        // Independent set variables
        public boolean inIndependentSet;
        public double independentSetWeight;
        public Interval independentSetPrevious;
        
        
        public Interval() {
            isSV=false; variantType=-1; chr=-1; first=-1; last=-1; weight=0.0; 
            vcfRecord=null; genotypes=null;
        }
        
        
        public final void clearIndependentSetVariables() {
            independentSetWeight=0.0; independentSetPrevious=null; inIndependentSet=false;
        }
        
        
        /**
         * @param record a VCF record;
         * @return FALSE iff $vcfRecord$ does not encode any known variant type.
         */
        public final boolean init(String record) {
            int i;
            int pos, length, nSamples;
            String tmpString;
            
            this.vcfRecord=record.split(VCF_SEPARATOR);
            tmpString=VCFconstants.getField(this.vcfRecord[7],VCFconstants.SVTYPE_STR);
            if (tmpString!=null && tmpString.length()!=0) {
                isSV=true;
                variantType=svType2Row(tmpString);
            }
            else {
                isSV=false;
                variantType=refAlt2Row(this.vcfRecord[3],this.vcfRecord[4]);
            }
            if (variantType==-1) return false;
            chr=VCFconstants.string2contig(this.vcfRecord[0]);
            pos=Integer.parseInt(this.vcfRecord[1]);
            tmpString=VCFconstants.getField(this.vcfRecord[7],VCFconstants.SVLEN_STR);
            if (tmpString!=null) {
                length=Integer.parseInt(tmpString);
                if (length<0) length=-length;
            }
            else if (variantType==REPLACEMENT) length=this.vcfRecord[3].length()-1;
            else length=Math.max(this.vcfRecord[3].length(),this.vcfRecord[4].length())-1;
            first=-1; last=-1;
            if (variantType==DEL || variantType==INV || variantType==DUP || variantType==REPLACEMENT) { first=pos+1; last=pos+length; }
            else if (variantType==INS) { first=pos; last=pos+1; }
            else if (variantType==SNP) { first=pos; last=pos; }
            nSamples=this.vcfRecord.length-9;
            if (genotypes==null || genotypes.length<nSamples) genotypes = new int[nSamples];
            for (i=0; i<nSamples; i++) genotypes[i]=gt2Id(this.vcfRecord[9+i]);
            return true;
        }
        
        
        /**
         * Sets $weight$ to a value loaded from the INFO field or the SAMPLE
         * field, depending on global variable $LOAD_WEIGHT_FROM_SAMPLE_COLUMN$.
         * If no value is found in the VCF record, $weight$ is set to one.
         */
        private final void setWeight(int sample) {
            int i, j, p;
            int gtLength;
            String gt, value;
            
            value=null;
            if (LOAD_WEIGHT_FROM_SAMPLE_COLUMN) {
                p=vcfRecord[8].indexOf(WEIGHT_ID_SAMPLE);
                if (p>=0) {
                    j=0;
                    for (i=0; i<p; i++) {
                        if (vcfRecord[8].charAt(i)==GT_SEPARATOR) j++;
                    }
                    gt=vcfRecord[9+sample]; gtLength=gt.length();
                    for (i=0; i<gtLength; i++) {
                        if (gt.charAt(i)!=GT_SEPARATOR) continue;
                        j--;
                        if (j>0) continue;
                        p=gt.indexOf(GT_SEPARATOR,i+1);
                        value=p>=0?gt.substring(i+1,p):gt.substring(i+1);
                        break;
                    }
                }
            }
            else value=VCFconstants.getField(vcfRecord[7],WEIGHT_ID_INFO);
            weight=value!=null?Double.parseDouble(value):1.0;
        }
        
        
        public final boolean isPresent(int sample) {
            return genotypes[sample]!=PHASED_00 && genotypes[sample]!=UNPHASED_00;
        }
        
        
        public final boolean onHap1(int sample) {
            return genotypes[sample]==PHASED_10 || genotypes[sample]==PHASED_11 || genotypes[sample]==UNPHASED_11;
        }
        
        
        public final boolean onHap2(int sample) {
            return genotypes[sample]==PHASED_01 || genotypes[sample]==PHASED_11 || genotypes[sample]==UNPHASED_11;
        }
        
        
        public final boolean isPhased(int sample) {
            return genotypes[sample]==PHASED_00 || genotypes[sample]==PHASED_01 || genotypes[sample]==PHASED_10 || genotypes[sample]==PHASED_11;
        }

        
        /**
         * @param nextInterval an interval that ends after $last$;
         * @param sample both this interval and $nextInterval$ are assumed to be
         * phased in $sample$;
         * @return TRUE iff this interval can precede $nextInterval$ in an 
         * independent set of $sample$.
         */
        public final boolean precedes(Interval nextInterval, int sample) {
            return N_GT_COLLISIONS[genotypes[sample]][nextInterval.genotypes[sample]]==0 || last<nextInterval.first;
        }
        
        
        public String toString() {
            return (isSV?"1":"0")+", "+variantType+", chr"+chr+"["+first+".."+last+"], weight="+weight;
        }
        
        
        /**
         * By last position only
         */
        public boolean equals(Object other) {
            final Interval otherInterval = (Interval)other;
            return last==otherInterval.last;
        }
        
        
        /**
         * By last position only
         */
        public int compareTo(Object other) {
            final Interval otherInterval = (Interval)other;
            if (last<otherInterval.last) return -1;
            else if (last>otherInterval.last) return 1;
            else return 0;
        }
        
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
    
    
    private static final void initType2color() {
        type2color = new Color[6];
        type2color[DEL] = new Color(0x006d9eeb);  // Blue
        type2color[INV] = new Color(0x0093c47d);  // Green
        type2color[DUP] = new Color(0x00ffd966);  // Yellow
        type2color[INS] = new Color(0x00cc0000);  // Red
        type2color[SNP] = new Color(0x00666666);  // Gray
        type2color[REPLACEMENT] = new Color(0x009900ff);  // Violet
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
                if (c=='.' || c=='0') return UNPHASED_00;
                else return UNPHASED_01;
            }
            else {
                if (c=='.' || c=='0') return UNPHASED_10;
                else return UNPHASED_11;
            }
        }
        else {
            if (a=='.' || a=='0') {
                if (c=='.' || c=='0') return PHASED_00;
                else return PHASED_01;
            }
            else {
                if (c=='.' || c=='0') return PHASED_10;
                else return PHASED_11;
            }
        }
    }
    
    
    
    
    
    
    
    
    
    
    


    // ---------- WRONG QUADRATIC IMPLEMENTATION OF INDEPENDENT SET. HAS COUNTEREXAMPLE. ------------
    // /**
//      * Remark: the procedure assumes that $window[0..last]$ is sorted by last
//      * position.
//      */
//     private static final void independentSet(int sample, int last) {
//         int i, j;
//         int sampleWindowLast;
//         double weight, currentWeight, maxWeight;
//         Interval currentInterval, neighbor;
//
//         // Computing the independent set
//         if (sampleWindow==null || sampleWindow.length<last+1) sampleWindow = new Interval[last+1];
//         sampleWindowLast=-1;
//         for (i=0; i<=last; i++) {
//             currentInterval=window[i];
//             if (!currentInterval.isPresent(sample) || !currentInterval.isPhased(sample)) continue;
//             sampleWindow[++sampleWindowLast]=currentInterval;
//             currentInterval.clearIndependentSetVariables();
//         }
//         maxWeight=0;
//         for (i=0; i<=sampleWindowLast; i++) {
//             currentInterval=sampleWindow[i];
//             currentWeight=currentInterval.weight;
//             for (j=i-1; j>=0; j--) {
//                 neighbor=sampleWindow[j];
//                 if (!neighbor.precedes(currentInterval,sample)) continue;
//                 weight=neighbor.independentSetWeight+currentWeight;
//                 if (weight>currentInterval.independentSetWeight) {
//                     currentInterval.independentSetWeight=weight;
//                     currentInterval.independentSetPrevious=neighbor;
// System.err.println("Setting independentSetWeight of "+currentInterval+" (GT="+currentInterval.genotypes[sample]+") to "+weight+" and independentSetPrevious="+neighbor+" (GT="+neighbor.genotypes[sample]+")");
//                 }
//             }
//             if (currentInterval.independentSetWeight>maxWeight) maxWeight=currentInterval.independentSetWeight;
//         }
//
//         // Removing from $sample$ every phased interval that is not in a
//         // max-weight independent set that ends farther to the right.
//         for (i=sampleWindowLast; i>=0; i--) {
//             if (sampleWindow[i].independentSetWeight!=maxWeight) continue;
//             currentInterval=sampleWindow[i];
//             while (currentInterval!=null) {
//                 currentInterval.inIndependentSet=true;
//                 currentInterval=currentInterval.independentSetPrevious;
//             }
//             break;
//         }
//         for (i=0; i<=sampleWindowLast; i++) {
//             if (!sampleWindow[i].inIndependentSet) sampleWindow[i].genotypes[sample]=UNPHASED_00;
//         }
//     }
    
}