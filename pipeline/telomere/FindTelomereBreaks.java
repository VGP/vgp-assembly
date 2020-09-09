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

public class FindTelomereBreaks {  
  private static final NumberFormat nf = new DecimalFormat("############.#");
  private static final int MIN_TEL = 24;

   public FindTelomereBreaks() {
   }

   public static void printUsage() {
      System.err.println("This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: getHist fasta1.fasta,fasta2.fasta");
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 3) { printUsage(); System.exit(1);}

      // initialize sizes
      HashMap<String, BitSet> scaffolds = new HashMap<String, BitSet>();
      HashMap<String, BitSet> finalScaffolds = new HashMap<String, BitSet>();
      HashMap<String, Integer> lengths = new HashMap<String, Integer>();
      BufferedReader bf = Utils.getFile(args[0], "lens");
      String line = null;
      while ((line = bf.readLine()) != null) {
         String[] split = line.trim().split("\\s+");
         if (!scaffolds.containsKey(split[0])) {
            scaffolds.put(split[0], new BitSet(Integer.parseInt(split[1])));
            finalScaffolds.put(split[0], new BitSet(Integer.parseInt(split[1])));
            lengths.put(split[0], Integer.parseInt(split[1]));
         }
      }
      bf.close();

      bf = Utils.getFile(args[1], "sdust");
      line = null;
      while ((line = bf.readLine()) != null) {
          String[] split = line.trim().split("\\s+");
          BitSet b = scaffolds.get(split[0]);
          b.set(Integer.parseInt(split[1]), Integer.parseInt(split[2]));
      }
      bf.close();

      // finally read the centromere and mark all those telomeric stretches in low complexity regions
      //
      // >asm_hic_scaffold_1	1083732140	0	46139	46145	6
      bf = Utils.getFile(args[2], "telomere");
      line = null;
      while ((line = bf.readLine()) != null) {
          String[] split = line.substring(1).trim().split("\\s+");
          BitSet b = scaffolds.get(split[0]);
          if (Integer.parseInt(split[5]) >= MIN_TEL) {
             int start = Math.max(0, Integer.parseInt(split[3])-100);
             int end = Math.min(lengths.get(split[0]), Integer.parseInt(split[4])+100);
             if (b.get(start, end).cardinality() >= (end - start)) {
                int rStart = (b.previousClearBit(start) < 0 ? 0 : b.previousClearBit(start));
                finalScaffolds.get(split[0]).set(rStart, b.nextClearBit(end));
             }
         }
      }
      bf.close();

      for (String scf : finalScaffolds.keySet()) {
         BitSet b = finalScaffolds.get(scf);

         for(int i = b.nextSetBit(0); i >= 0; i = b.nextSetBit(i+1)) {
            int end = b.nextClearBit(i) - 1;
            System.out.println("Found telomere positiosn " + i + " to " + end + " is a telomere in " + scf + " of length " + lengths.get(scf));
            i = end;
         }
      }
   }
}
