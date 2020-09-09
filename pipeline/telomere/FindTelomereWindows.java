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
      System.err.println("Usage: java -jar FindTelomereWindows.jar <in> <identity> <threshold>");
      System.err.println("This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: getHist fasta1.fasta,fasta2.fasta");
      System.err.println("");
   }
   public static void processScaffold(String name, BitSet b, int length) {
      if (b == null) { return; }


       for (int i = MIN_OFFSET; i <= length; i+=WINDOW_SIZE/5) {
         int car = b.get(i, Math.min(length, i+WINDOW_SIZE)).cardinality();
         int den = Math.min(WINDOW_SIZE, length-i);
         if ((double)car / den >= THRESHOLD)
            System.out.println("Window\t" + name + "\t" + length + "\t" + i + "\t" + (i+den) + "\t" + ((double)car / den));

         if (i+WINDOW_SIZE >= length)
            break;
     }
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}

      if (args.length == 3) {
          THRESHOLD = Double.parseDouble(args[2]);
      }

      // initialize sizes
      BitSet scaffold = null;
      String name = null;
      int length = 0;
      Double identity = Double.parseDouble(args[1]) / 100;
      THRESHOLD = THRESHOLD * Math.pow(identity, 6);
      System.err.println("Given error rate of " + identity + " running with adjusted threshold of " + THRESHOLD + " due to survival prob " + Math.pow(identity, 6));

      BufferedReader bf = Utils.getFile(args[0], "telomere");
      String line = null;
      while ((line = bf.readLine()) != null) {
          String[] split = line.trim().split("\\s+");
          if (scaffold == null || !split[0].equalsIgnoreCase(name)) {
             processScaffold(name, scaffold, length);
             length = Integer.parseInt(split[split.length-5]);
             scaffold = new BitSet(length);
             name = split[0];
          }
          //ignoring strandedness for now
          scaffold.set(Integer.parseInt(split[split.length-3]), Integer.parseInt(split[split.length-2]));
      }
      processScaffold(name, scaffold, length);
      bf.close();
   }
}
