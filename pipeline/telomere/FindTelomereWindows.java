import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.io.File;
import java.io.PrintStream;
import java.util.BitSet;

public class FindTelomereWindows {  
  private static final NumberFormat nf = new DecimalFormat("############.#");
  private static final int WINDOW_SIZE = 1000;
  private static final int MIN_OFFSET = 0;
  private static double THRESHOLD = 0.4;

   public FindTelomereWindows() {
   }

   public static void printUsage() {
      System.err.println("This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: getHist fasta1.fasta,fasta2.fasta");
   }
   public static void processScaffold(String name, BitSet b) {
      if (b == null) { return; }
      for (int i = MIN_OFFSET; i <= b.size()-MIN_OFFSET-WINDOW_SIZE; i+=WINDOW_SIZE/5) {
        int car = b.get(i, i+WINDOW_SIZE).cardinality();
        if ((double)car / WINDOW_SIZE >= THRESHOLD) {
           System.out.println("Window\t" + name + "\t" + b.size() + "\t" + i + "\t" + (i+WINDOW_SIZE) + "\t" + ((double)car / WINDOW_SIZE));
        }
     }
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}

      // initialize sizes
      BitSet scaffold = null;
      String name = null;
      Double identity = Double.parseDouble(args[1]) / 100;
      THRESHOLD = THRESHOLD * Math.pow(identity, 6);
      System.err.println("Given error rate of " + identity + " running with adjusted threshold of " + THRESHOLD + " due to survival prob " + Math.pow(identity, 6));

      BufferedReader bf = Utils.getFile(args[0], "telomere");
      String line = null;
      while ((line = bf.readLine()) != null) {
          String[] split = line.trim().split("\\s+");
          if (scaffold == null || !split[0].equalsIgnoreCase(name)) {
             processScaffold(name, scaffold);
             scaffold = new BitSet(Integer.parseInt(split[1]));
             name = split[0];
          }
          //ignoring strandedness for now
          scaffold.set(Integer.parseInt(split[3]), Integer.parseInt(split[4]));
      }
      processScaffold(name, scaffold);
      bf.close();
   }
}
