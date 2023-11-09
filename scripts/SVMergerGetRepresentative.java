import java.util.Arrays;
import java.util.HashMap;
import java.io.*;


/**
 * 
 */
public class SVMergerGetRepresentative {
    /**
     * @param args 0: assumed to be sorted by clique ID.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        
        int id;
        String str, currentClique;
        BufferedReader br;
        String[] tokens;
        HashMap<Pair,Pair> map;
        Pair key, value, newPair;
        Pair[] pairs;
        
        key = new Pair();
        map = new HashMap<Pair,Pair>();
        br = new BufferedReader(new FileReader(INPUT_FILE));
        currentClique="";
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            if (currentClique.length()==0) currentClique=tokens[8];
            if (tokens[8].equals(currentClique)) {
                key.start=Integer.parseInt(tokens[2]);
                key.end=Integer.parseInt(tokens[3]);
                value=map.get(key);
                if (value==null) {
                    newPair = new Pair(tokens[0],tokens[1],key.start,key.end);
                    map.put(newPair,newPair);
                }
                else value.count++;
            }
            else {
                pairs = new Pair[map.size()];
                map.keySet().toArray(pairs);
                if (pairs.length>1) Arrays.sort(pairs);
                System.out.println(pairs[0].id);
                currentClique=tokens[8];
                map.clear();
                newPair = new Pair(tokens[0],tokens[1],Integer.parseInt(tokens[2]),Integer.parseInt(tokens[3]));
                map.put(newPair,newPair);
            }
            str=br.readLine();
        }
        pairs = new Pair[map.size()];
        map.keySet().toArray(pairs);
        if (pairs.length>1) Arrays.sort(pairs);
        System.out.println(pairs[0].id);
        br.close();
    }
    
    
    public static class Pair implements Comparable {
        public String id;
        public String chr;
        public int start, end;
        public int count;
        
        public Pair() { }
        
        public Pair(String id, String chr, int start, int end) {
            this.id=id; this.chr=chr;
            this.start=start; this.end=end;
            count=1;
        }
        
        public int compareTo(Object other) {
            Pair otherPair = (Pair)other;
            if (count>otherPair.count) return -1;
            else if (count<otherPair.count) return 1;
            return 0;
        }
        
        public boolean equals(Object other) {
            Pair otherPair = (Pair)other;
            return otherPair.start==start && otherPair.end==end;
        }
        
        public int hashCode() {
            return (start+","+end).hashCode();
        }
    }
    
}