import java.io.*;


public class OptimizerTimeToSeconds {

    public static void main(String[] args) throws Exception {
        final String INPUT_FILE = args[0];
        
        int p;
        double value;
        String str;
        BufferedReader br;
        
        br = new BufferedReader(new FileReader(INPUT_FILE));
        str=br.readLine();
        while (str!=null) {
            value=0.0;
            p=str.indexOf("h");
            if (p>=0) {
                value+=Double.parseDouble(str.substring(0,p))*60*60;
                str=str.substring(p+1);
            }
            p=str.indexOf("m");
            if (p>=0 && p+1<str.length() && str.charAt(p+1)!='s') {
                value+=Double.parseDouble(str.substring(0,p))*60;
                str=str.substring(p+1);
            }
            p=str.indexOf("s");
            if (p>=0 && p>0 && str.charAt(p-1)!='m' && str.charAt(p-1)!='u') {
                value+=Double.parseDouble(str.substring(0,p));
                str=str.substring(p+1);
            }
            p=str.indexOf("ms");
            if (p<0) p=str.indexOf("us");
            if (p>=0) value+=Double.parseDouble(str.substring(0,p))/1000;
            System.out.println(((int)value)+"");
            str=br.readLine();
        }
        br.close();
    }

}