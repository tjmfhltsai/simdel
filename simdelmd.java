import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import java.io.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.image.*;
import java.util.*;
class DNA
{
    String id;
    String sec;
}
public class simdelmd 
{
    static double ratio=0.8;
    int filept;
    static String fname="/mnt/dos/bio/su/aa",sfname="/mnt/dos/bio/su/bb";
    static String outdir="/mnt/dos/bio/su/out";
    String logname="log.txt",sourcefile;
    static int seglen=100;
    String outfname;
    Vector<DNA> dnalist=new Vector<DNA>();
    Thread t;
    public simdelmd(String source, String outf,String log)   
    {
        
        logname=log;
        outfname=outf;
        sourcefile=source;
        t= new Thread() {
                  public void run() {
                        go();
               }
        };
        t.start();
    
    }

    public void go(){
        clearlog();
        if ((new File(sourcefile)).isFile()) {
            readdna(sourcefile);
            log("read "+sourcefile);
            File file = new File(fname);
            if (file.isDirectory()){ //each input secquence
                File[] listOfFiles = file.listFiles();
                for (File filei : listOfFiles) {
                    if (filei.isFile()) {   
                        log("compare with "+filei.getPath()+"  similarity="+ratio);
                        proc(filei.getPath());
                    }
                }
            }
                log("write output to "+outfname);
                outfile(outfname);
        }
          
        
    }
    void clearlog()
   {
       try{ 
           BufferedWriter sbw=new BufferedWriter(new FileWriter(logname));
           sbw.write("");         
           sbw.flush();
           sbw.close();
       }
        catch (Exception ee){System.out.println(ee);}
   }
    void log(String msg)
   {
       try{ 
           BufferedWriter sbw=new BufferedWriter(new FileWriter(logname,true));
           sbw.write(msg+"\n");         
           sbw.flush();
           sbw.close();
       }
        catch (Exception ee){System.out.println(ee);}
   }
   void outfile(String sfn)
   {
       try{
           
           BufferedWriter sbw=new BufferedWriter(new FileWriter(sfn));
           String tmp,sstr="";
           DNA now;
           for (int i=0;i<dnalist.size();i++){
               now=dnalist.elementAt(i);
               sbw.write(now.id+"\n");
               int pt=0,pte=60;
               while (pt<now.sec.length())
               {
                   if (pt+60>now.sec.length()) pte=now.sec.length();
                   else pte=pt+60;
                   sbw.write(now.sec.substring(pt,pte)+"\n");
                   pt=pte;
               }
           }
           sbw.flush();
           sbw.close();
       }
        catch (Exception ee){System.out.println(ee);}
   }
   void readdna(String sfn)
   {
       try{
           DNA now=new DNA();
           BufferedReader sbr=new BufferedReader(new FileReader(sfn));
           String tmp,sstr="";
           boolean first=true;
           while ((tmp=sbr.readLine())!=null)
           {
               tmp=tmp.trim();
               if (tmp.length()>=6 && tmp.substring(0,6).equals(">NODE_"))
               {
                   if (sstr.length()>0){
                       now.sec=sstr;
                       //System.out.println(now.id);
                       dnalist.add(now);
                       now=new DNA();
                       now.id=tmp;
                       sstr="";
                   }else{
                       if (first){
                           now.id=tmp;
                           first=false;
                       }
                   }
                   continue;
               }
               sstr=sstr+tmp;
           }
           if (sstr.length()>0){
                       now.sec=sstr;
                       //System.out.println(now.id);
                       dnalist.add(now);
                       
           }
           sbr.close();
       }
        catch (Exception ee){System.out.println(ee);}
   }
   void proc(String fn)
   {
   
       try{
           System.out.println(fn);    
           BufferedReader br=new BufferedReader(new FileReader(fn)); 
           String str="",tmp;
           boolean beg=false;
           while ((tmp=br.readLine())!=null)
           {
               if (!beg && tmp.indexOf("\\par")>=0)
               {
                   if (!beg) {
                       beg=true;
                       continue;
                   }    
               }
               tmp=tmp.trim();
               if (beg && tmp.length()==0){
                   beg=false;
                   break;
               }
               if (beg)  str=str+tmp.replaceAll("\\\\par","");
           }
           
           br.close();
           //System.out.println(dnalist.size());
           for (int i=0;i<dnalist.size();i++)
           {
               if (ratio==1)
                   match(dnalist.elementAt(i),str);
               else {
                   pmatch(dnalist.elementAt(i),str);
               }
           }
       }
       catch (Exception ee){System.out.println(ee);}
   }
   void match(DNA dna,String str)
   {
       int pt;
       pt=dna.sec.indexOf(str);
       if (pt>=0) {
           System.out.println("hit");
           log(dna.id+" hit position:"+pt+"  "+str);
           dna.sec=dna.sec.replaceAll(str,"");
       }
       
   }
   void pmatch(DNA dna,String str)
   {
       //System.out.println("aaa="+dna.id);
       if (dna.sec.length()<str.length()) return;
       int pt=0,range=(int)((1-ratio)*str.length());//search area
       int i,last=dna.sec.length()-str.length()+1;
       int len=str.length();
       double sim=0,max=0;
       int maxi=-1;
       try{
       for (i=0;i<last;i++)
       {
           //System.out.println("i="+i+" last="+last);
           sim=similarity(dna.sec.substring(i,i+len), str);
           if (sim>ratio)//hit
           {
               max=sim;
               maxi=i;
               pt=1;
               while (i+pt<last && pt<=range)//search max similar
               {
                   sim=similarity(dna.sec.substring(i+pt,i+pt+len), str);
                   if (sim>max){
                       max=sim;
                       maxi=i+pt;
                   }
                   pt++;
               }
               maxi=maxi+1;
               log(dna.id+" hit position:"+maxi+"  similarity:"+max+" pattern:"+str);
               //log("befor="+dna.sec);
               if (maxi>=0)
		{
			if (dna.sec.length()>maxi+len)
			{
               			dna.sec=dna.sec.substring(0,maxi)+
					dna.sec.substring(maxi+len);
			}
			else
               			dna.sec=dna.sec.substring(0,maxi);
               		last=dna.sec.length()-str.length()+1;
               		i=maxi-2;
		}
           }
           
       }
      }
      catch (Exception ex){
	System.out.println("maxi="+maxi+" len="+len);
        ex.printStackTrace();
      }
     }
  /**
   * Calculates the similarity (a number within 0 and 1) between two strings.
   */
  public static double similarity(String s1, String s2) {
    String longer = s1, shorter = s2;
    if (s1.length() < s2.length()) { // longer should always have greater length
      longer = s2; shorter = s1;
    }
    int longerLength = longer.length();
    if (longerLength == 0) { return 1.0; /* both strings are zero length */ }
    /* // If you have Apache Commons Text, you can use it to calculate the edit distance:
    LevenshteinDistance levenshteinDistance = new LevenshteinDistance();
    return (longerLength - levenshteinDistance.apply(longer, shorter)) / (double) longerLength; */
    int stopcnt=(int)(longerLength*(1-ratio));
    int total=0;
    while (shorter.length()>=2*seglen)
    {
        String test=shorter.substring(0,seglen);
        shorter=shorter.substring(seglen);
        String test1=longer.substring(0,seglen);
        longer=longer.substring(seglen);
    	int stopcnt1=(int)(seglen*(1-ratio));
        int dist=getLevenshteinDistance(test1, test,stopcnt1*3/2);
        if (dist<0)
            total=total+test.length();
        else total=total+dist;
        //System.out.println("totla="+total+" "+shorter.length());
        if (total>stopcnt) return 0;
        
    }
    if (shorter.length()>0)
    {
    	int stopcnt1=(int)(shorter.length()*(1-ratio));
        int dist=getLevenshteinDistance(longer, shorter,stopcnt1*3/2);
        if (dist<0)
            total=total+shorter.length();
        else total=total+dist;
        if (total>stopcnt) return 0;
    }
    
        return (longerLength - total) / (double) longerLength;
    
  }

  // Example implementation of the Levenshtein Edit Distance
  // See http://rosettacode.org/wiki/Levenshtein_distance#Java
  public static int editDistance(String s1, String s2) {
    s1 = s1.toLowerCase();
    s2 = s2.toLowerCase();

    int[] costs = new int[s2.length() + 1];
    for (int i = 0; i <= s1.length(); i++) {
      int lastValue = i;
      for (int j = 0; j <= s2.length(); j++) {
        if (i == 0)
          costs[j] = j;
        else {
          if (j > 0) {
            int newValue = costs[j - 1];
            if (s1.charAt(i - 1) != s2.charAt(j - 1))
              newValue = Math.min(Math.min(newValue, lastValue),
                  costs[j]) + 1;
            costs[j - 1] = lastValue;
            lastValue = newValue;
          }
        }
      }
      if (i > 0)
        costs[s2.length()] = lastValue;
    }
    return costs[s2.length()];
  }
  /**
 
 * StringUtils.getLevenshteinDistance("","", 0)               = 0
 * StringUtils.getLevenshteinDistance("aaapppp", "", 8)       = 7
 * StringUtils.getLevenshteinDistance("aaapppp", "", 7)       = 7
 * StringUtils.getLevenshteinDistance("aaapppp", "", 6))      = -1
 * StringUtils.getLevenshteinDistance("elephant", "hippo", 7) = 7
 * StringUtils.getLevenshteinDistance("elephant", "hippo", 6) = -1
 * StringUtils.getLevenshteinDistance("hippo", "elephant", 7) = 7
 * StringUtils.getLevenshteinDistance("hippo", "elephant", 6) = -1
 * </pre>
 *
 * @param s  the first String, must not be null
 * @param t  the second String, must not be null
 * @param threshold the target threshold, must not be negative
 * @return result distance, or {@code -1} if the distance would be greater than the threshold
 * @throws IllegalArgumentException if either String input {@code null} or negative threshold
 */
public static int getLevenshteinDistance(CharSequence s, CharSequence t, final int threshold) {
    if (s == null || t == null) {
        throw new IllegalArgumentException("Strings must not be null");
    }
    if (threshold < 0) {
        throw new IllegalArgumentException("Threshold must not be negative");
    }

    /*
    This implementation only computes the distance if it's less than or equal to the
    threshold value, returning -1 if it's greater.  The advantage is performance: unbounded
    distance is O(nm), but a bound of k allows us to reduce it to O(km) time by only
    computing a diagonal stripe of width 2k + 1 of the cost table.
    It is also possible to use this to compute the unbounded Levenshtein distance by starting
    the threshold at 1 and doubling each time until the distance is found; this is O(dm), where
    d is the distance.

    One subtlety comes from needing to ignore entries on the border of our stripe
    eg.
    p[] = |#|#|#|*
    d[] =  *|#|#|#|
    We must ignore the entry to the left of the leftmost member
    We must ignore the entry above the rightmost member

    Another subtlety comes from our stripe running off the matrix if the strings aren't
    of the same size.  Since string s is always swapped to be the shorter of the two,
    the stripe will always run off to the upper right instead of the lower left of the matrix.

    As a concrete example, suppose s is of length 5, t is of length 7, and our threshold is 1.
    In this case we're going to walk a stripe of length 3.  The matrix would look like so:

       1 2 3 4 5
    1 |#|#| | | |
    2 |#|#|#| | |
    3 | |#|#|#| |
    4 | | |#|#|#|
    5 | | | |#|#|
    6 | | | | |#|
    7 | | | | | |

    Note how the stripe leads off the table as there is no possible way to turn a string of length 5
    into one of length 7 in edit distance of 1.

    Additionally, this implementation decreases memory usage by using two
    single-dimensional arrays and swapping them back and forth instead of allocating
    an entire n by m matrix.  This requires a few minor changes, such as immediately returning
    when it's detected that the stripe has run off the matrix and initially filling the arrays with
    large values so that entries we don't compute are ignored.

    See Algorithms on Strings, Trees and Sequences by Dan Gusfield for some discussion.
     */

    int n = s.length(); // length of s
    int m = t.length(); // length of t

    // if one string is empty, the edit distance is necessarily the length of the other
    if (n == 0) {
        return m <= threshold ? m : -1;
    } else if (m == 0) {
        return n <= threshold ? n : -1;
    }

    if (n > m) {
        // swap the two strings to consume less memory
        final CharSequence tmp = s;
        s = t;
        t = tmp;
        n = m;
        m = t.length();
    }

    int p[] = new int[n + 1]; // 'previous' cost array, horizontally
    int d[] = new int[n + 1]; // cost array, horizontally
    int _d[]; // placeholder to assist in swapping p and d

    // fill in starting table values
    final int boundary = Math.min(n, threshold) + 1;
    for (int i = 0; i < boundary; i++) {
        p[i] = i;
    }
    // these fills ensure that the value above the rightmost entry of our
    // stripe will be ignored in following loop iterations
    Arrays.fill(p, boundary, p.length, Integer.MAX_VALUE);
    Arrays.fill(d, Integer.MAX_VALUE);

    // iterates through t
    for (int j = 1; j <= m; j++) {
        final char t_j = t.charAt(j - 1); // jth character of t
        d[0] = j;

        // compute stripe indices, constrain to array size
        final int min = Math.max(1, j - threshold);
        final int max = (j > Integer.MAX_VALUE - threshold) ? n : Math.min(n, j + threshold);

        // the stripe may lead off of the table if s and t are of different sizes
        if (min > max) {
            return -1;
        }

        // ignore entry left of leftmost
        if (min > 1) {
            d[min - 1] = Integer.MAX_VALUE;
        }

        // iterates through [min, max] in s
        for (int i = min; i <= max; i++) {
            if (s.charAt(i - 1) == t_j) {
                // diagonally left and up
                d[i] = p[i - 1];
            } else {
                // 1 + minimum of cell to the left, to the top, diagonally left and up
                d[i] = 1 + Math.min(Math.min(d[i - 1], p[i]), p[i - 1]);
            }
        }

        // copy current distance counts to 'previous row' distance counts
        _d = p;
        p = d;
        d = _d;
    }

    // if p[n] is greater than the threshold, there's no guarantee on it being the correct
    // distance
    if (p[n] <= threshold) {
        return p[n];
    }
    return -1;
}
    public static void main(String args[]) 
    {
        if (args.length<5)
        {
            System.out.println("java simdel similarity(0.8-1.0) input_dir database_dir output_dir segment_length");
            System.exit(0);
        }
        ratio=Double.parseDouble(args[0]);
        fname=args[1];
        sfname=args[2];
        outdir=args[3];
        seglen=Integer.parseInt(args[4]);
        File sfile = new File(sfname);  
        File[] listOfsFiles = sfile.listFiles();
        simdelmd app[]=new simdelmd[listOfsFiles.length];
        int pp=0;
            for (File sfilei : listOfsFiles) { //for each database secquence
                if (sfilei.isFile()) {
                    String outfname=outdir+"/"+sfilei.getName();
                    String logname=outdir+"/"+sfilei.getName()+".log";
                    app[pp]=new simdelmd(sfilei.getPath(),outfname,logname);
                    pp++;
                    }
                    
                }//if (sfilei.isFile())       
    }
}
