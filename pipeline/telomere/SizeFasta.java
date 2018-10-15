import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class SizeFasta {  
  private static final NumberFormat nf = new DecimalFormat("############.#");
   
   private static final String[] FASTA_ENDS = {"trimmed", "contig", "RELEASE5", "bases", "qual", "fasta", "fas", "fna", "fa", "genome"};
   private static final String[] FASTQ_ENDS = {"seq", "txt", "fastq", "fq"};
   private static final int MAX_DEFAULT_STRING = 4000000;
   private boolean skipNs = false;
   private boolean stopEarly = false;

   public SizeFasta() {
   }

   public void processFasta(String inputFile) throws Exception {
      BufferedReader bf = Utils.getFile(inputFile, FASTA_ENDS, true);
      
      String line = null;
      StringBuffer fastaSeq = new StringBuffer(MAX_DEFAULT_STRING);
      String header = "";
      
      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) {
            if (header.length() > 0) {
System.out.println(header + "\t" + (skipNs ? fastaSeq.toString().replaceAll("N", "").length() : fastaSeq.length()));
if (stopEarly) { System.exit(0); }
            }
            header = line.substring(1);
            fastaSeq.setLength(0);
            //fastaSeq = new StringBuffer();
         }
         else {
            fastaSeq.append(line);
/*
for (int i = 0; i < line.length(); i++) {
if (line.charAt(i) != 'A' && line.charAt(i) != 'C' && line.charAt(i) != 'G' && line.charAt(i) != 'T') {
System.err.println("Error invalid line " + line + " in header " + header);
System.exit(1);
}
}
*/
         }
      }

      if (fastaSeq.length() != 0) { System.out.println(header + "\t" + (skipNs ? fastaSeq.toString().replaceAll("N", "").length() : fastaSeq.length()));}
      bf.close();
   }

   public long sizeSingleRecord(String inputFile) throws Exception {
      BufferedReader bf = Utils.getFile(inputFile, FASTA_ENDS);

      String line = null;
      StringBuffer fastaSeq = new StringBuffer(MAX_DEFAULT_STRING);
      String header = "";
      long totalLen = 0;

      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) {
            if (header.length() > 0) totalLen+=fastaSeq.length();
            header = line.substring(1);
            fastaSeq.setLength(0); // =new StringBuffer();
         }
         else {
            fastaSeq.append(line);
         }
      }

      if (fastaSeq.length() != 0) { totalLen+=fastaSeq.length();}

      bf.close();
      return totalLen;
   }

   public void processFastq(String inputFile) throws Exception {
      BufferedReader bf = Utils.getFile(inputFile, FASTQ_ENDS);
if (bf == null) {
bf = Utils.getFile(inputFile, "fq");
}
      String line = null;

      while ((line = bf.readLine()) != null) {
         // read four lines at a time for fasta, qual, and headers
         String ID = line.split("\\s+")[0].substring(1).trim();
         String fasta = bf.readLine();
String qualLine = bf.readLine();
         String qualID = qualLine.split("\\s+")[0].substring(1).trim();

         if (qualID.length() != 0 && !qualID.equals(ID)) {
            System.err.println("Error ID " + ID + " (" + qualID.length() + ") DOES not match quality ID " + qualID);
            System.exit(1);
         }
         String qualSeq = bf.readLine();
         System.out.println(ID + "\t" + qualSeq.length() + "\t" + (skipNs ? fasta.replaceAll("N", "").length() : fasta.length()));
         if (stopEarly) { System.exit(0); }
      }

      bf.close();
   }

   public static void printUsage() {
      System.err.println("This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: SizeFasta fasta1.fasta,fasta2.fasta");
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}

      SizeFasta f = new SizeFasta();

      for (int i = 0; i < args.length; i++) {
        if (args[i].equalsIgnoreCase("-skip")) { f.skipNs = true; continue; }
        if (args[i].equalsIgnoreCase("-first")) {f.stopEarly = true; continue; }
        String[] splitLine = args[i].trim().split(",");
        for (int j = 0; j < splitLine.length; j++) {
          boolean done = false;
System.err.println("Processing file " + splitLine[j]);
          for (String s : FASTA_ENDS) {
     	    if (splitLine[j].contains(s) && !splitLine[j].contains("fastq") && !splitLine[j].contains("fq") && !splitLine[j].contains("txt")) {
               f.processFasta(splitLine[j]);
               done = true;
               break;
            }
          }
          if (!done) {
             if (splitLine[j].contains("seq") || splitLine[j].contains("fastq") || splitLine[j].contains("fq") | splitLine[j].contains("txt")) {
                f.processFastq(splitLine[j]);
             } else {
                // assume fasta if not known
                f.processFasta(splitLine[j]);
                //System.err.println("Unknown file type " + splitLine[j]);
             }
          }
       }
     }
   }
}
